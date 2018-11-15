#include "udf.h"
#include "sg_mphase.h"
#include "mem.h"
#include "metric.h"
#include "flow.h"
#include "sg.h"

#include "thermoProp.h"
#include "funDecs.h"

#define TLAMBDA 2.1768
#define TREF 1.4
#define MAG_GRAD_TP_CONV_CRIT 1e-2
#define MAG_GRAD_TP_VALUE 1e-2

/*
   Description of the user defined memories:

   CUDMI(c,t,0) -> x component of the vector (-gradP/rho/s + gradT) 
   CUDMI(c,t,1) -> y component of the vector (-gradP/rho/s + gradT) 
   CUDMI(c,t,2) -> z component of the vector (-gradP/rho/s + gradT) 
   CUDMI(c,t,3) -> |G| = magnitude of the vector (-gradP/rho/s + gradT) 
   CUDMI(c,t,4) -> thermal conductivity 
   CUDMI(c,t,5) -> density
   CUDMI(c,t,6) -> entropy
   CUDMI(c,t,7) -> dynamic viscosity
   CUDMI(c,t,8) -> Gorter-Mellink parameter
   CUDMI(c,t,9) -> rho_n
   CUDMI(c,t,10) -> rho_s
   CUDMI(c,t,11) -> B = (s/A/rho_n/magGradTp^2)^(1/3)
   CUDMI(c,t,12) -> S1x = divV_x
   CUDMI(c,t,13) -> S1y = divV_y
   CUDMI(c,t,14) -> S1z = divV_z
   CUDMI(c,t,15) -> x component of laplacian(-gradp/rho/s + gradT)
   CUDMI(c,t,16) -> y component of laplacian(-gradp/rho/s + gradT)
   CUDMI(c,t,17) -> z component of laplacian(-gradp/rho/s + gradT)
   CUDMI(c,t,18) -> x component of 1/3grad(div(-gradp/rho/s + gradT))
   CUDMI(c,t,19) -> y component of 1/3grad(div(-gradp/rho/s + gradT))
   CUDMI(c,t,20) -> z component of 1/3grad(div(-gradp/rho/s + gradT))
   CUDMI(c,t,21) -> z component of boundary condition Eq. (15)

   FUDMI(f,t,22) -> rho_s/rho*(s/A/rho_n/magGradTp^2)^(1/3)*gradTp_z 
*/

enum
{
	Gx1,      /* x component of the vector (-1/s/rho*gradP + gradT) */
	Gx2,      /* y component of the vector (-1/s/rho*gradP + gradT) */
	Gx3,      /* z component of the vector (-1/s/rho*gradP + gradT) */
	GGx11,    /* element x11 of diad GG */
	GGx12,    /* element x12 of diad GG */
	GGx13,    /* element x13 of diad GG */
	GGx21,    /* element x21 of diad GG */
	GGx22,    /* element x22 of diad GG */
	GGx23,    /* element x23 of diad GG */
	GGx31,    /* element x31 of diad GG */
	GGx32,    /* element x32 of diad GG */
	GGx33,    /* element x33 of diad GG */
	gradGx11, /* x11 component of the tensor T2 = grad(G) = grad(-1/s/rho*gradP + gradT) */
	gradGx12, /* x12 component of the tensor T2 = grad(G) = grad(-1/s/rho*gradP + gradT) */
	gradGx13, /* x13 component of the tensor T2 = grad(G) = grad(-1/s/rho*gradP + gradT) */
	gradGx21, /* x21 component of the tensor T2 = grad(G) = grad(-1/s/rho*gradP + gradT) */
	gradGx22, /* x22 component of the tensor T2 = grad(G) = grad(-1/s/rho*gradP + gradT) */
	gradGx23, /* x23 component of the tensor T2 = grad(G) = grad(-1/s/rho*gradP + gradT) */
	gradGx31, /* x31 component of the tensor T2 = grad(G) = grad(-1/s/rho*gradP + gradT) */
	gradGx32, /* x32 component of the tensor T2 = grad(G) = grad(-1/s/rho*gradP + gradT) */
	gradGx33, /* x33 component of the tensor T2 = grad(G) = grad(-1/s/rho*gradP + gradT) */
	divG,     /* divergence of the G vector */

	v0e,  /* k*gradTp_x in energy equation */
	v1e,  /* k*gradTp_y in energy equation */
	v2e,  /* k*gradTp_z in energy equation */
	vx1x2, /* x1x2 component of v dyad */
	vx1y2, /* x1y2 component of v dyad */
	vx1z2, /* x1z2 component of v dyad */
	vy1x2, /* y1x2 component of v dyad */
	vy1y2, /* y1y2 component of v dyad */
	vy1z2, /* y1z2 component of v dyad */
	vz1x2, /* z1x2 component of v dyad */
	vz1y2, /* z1y2 component of v dyad */
	vz1z2, /* z1z2 component of v dyad */
	N_REQUIRED_UDS
};

/* declaration of the 3D vector */
struct vector3d
{
	double x1, x2, x3;
};

/* calculates the magnitude of a 3D vector */
double magVector3d(struct vector3d v)
{
	return sqrt(v.x1*v.x1 + v.x2*v.x2 + v.x3*v.x3);
}

/* declaration of the dyad of the second rank */
struct dyad2rank
{
	double x11, x12, x13,
	     x21, x22, x23,
		 x31, x32, x33;
};

/* calculates the dyadic product of two vectors */
struct dyad2rank dyadicProd(const struct vector3d u, const struct vector3d v)
{
	struct dyad2rank d;
	/* the first row */
	d.x11 = u.x1*v.x1;
	d.x12 = u.x1*v.x2;
	d.x13 = u.x1*v.x3;

	/* the second row */
	d.x21 = u.x2*v.x1;
	d.x22 = u.x2*v.x2;
	d.x23 = u.x2*v.x3;

	/* the third row */
	d.x31 = u.x3*v.x1;
	d.x32 = u.x3*v.x2;
	d.x33 = u.x3*v.x3;

	return d;
}

/* calculates product of a dyad and a constant */
struct dyad2rank dyadMultByConst(struct dyad2rank d, const double con)
{
	d.x11 *= con;
	d.x12 *= con;
	d.x13 *= con;

	d.x21 *= con;
	d.x22 *= con;
	d.x23 *= con;

	d.x31 *= con;
	d.x32 *= con;
	d.x33 *= con;

	return d;
}

/* calculates the vector: (-gradP/rho/s + gradT) = G */
/* it runs at the begining of the ADJUST macro */
struct vector3d gradTp(cell_t c, Thread *t)
{
	struct vector3d gradp;
	struct vector3d gradT;
	double magSumGradpbyrhosGradT;

	C_UDMI(c,t,6) = sArray(c,t);
	double s = C_UDMI(c,t,6);
	double rho = C_UDMI(c,t,5); /* already calculated in thermophysical properties */

	gradp.x1 = C_P_G(c,t)[0];
	gradp.x2 = C_P_G(c,t)[1];
	gradp.x3 = C_P_G(c,t)[2];

	gradT.x1 = C_T_G(c,t)[0];
	gradT.x2 = C_T_G(c,t)[1];
	gradT.x3 = C_T_G(c,t)[2];
	
	struct vector3d sumGradpbyrhosGradT;
	sumGradpbyrhosGradT.x1 = -gradp.x1/rho/s + gradT.x1;
	sumGradpbyrhosGradT.x2 = -gradp.x2/rho/s + gradT.x2;
	sumGradpbyrhosGradT.x3 = -gradp.x3/rho/s + gradT.x3;

	C_UDSI(c,t,Gx1) = sumGradpbyrhosGradT.x1;
	C_UDSI(c,t,Gx2) = sumGradpbyrhosGradT.x2;
	C_UDSI(c,t,Gx3) = sumGradpbyrhosGradT.x3;

	/* filling the components of the tensor T2 = grad(G) = grad(-1/rho/s*gradp + gradT) */
	/* the first row */
	C_UDSI(c, t, gradGx11) = C_UDSI_G(c, t, Gx1)[0];
	C_UDSI(c, t, gradGx12) = C_UDSI_G(c, t, Gx2)[0];
	C_UDSI(c, t, gradGx13) = C_UDSI_G(c, t, Gx3)[0];

	/* the second row */
	C_UDSI(c, t, gradGx21) = C_UDSI_G(c, t, Gx1)[1];
	C_UDSI(c, t, gradGx22) = C_UDSI_G(c, t, Gx2)[1];
	C_UDSI(c, t, gradGx23) = C_UDSI_G(c, t, Gx3)[1];

	/* the third row */
	C_UDSI(c, t, gradGx31) = C_UDSI_G(c, t, Gx1)[2];
	C_UDSI(c, t, gradGx32) = C_UDSI_G(c, t, Gx2)[2];
	C_UDSI(c, t, gradGx33) = C_UDSI_G(c, t, Gx3)[2];

	return sumGradpbyrhosGradT;
}

/* calculates the div of v dyad */
void divDyadV(cell_t c, Thread *t)
{
	/* rhon, rhos and B = C_UDMI(c,t,11) are already calculated in ADJUST calcSrcTermsCmpts */
	double rhon = C_UDMI(c,t,9);
	double rhos = C_UDMI(c,t,10);
	/* rho is already calculated during calculating properties by the solver */
	/* it is done at the end of each iteration */
	double rho = C_UDMI(c,t,5);
	double con = rhon*rhos*C_UDMI(c,t,11)*C_UDMI(c,t,11)/rho;

	/* calculated at the begining of the ADJUST calcSrcTermsCmpts */
	struct vector3d GradTp;
	GradTp.x1 = C_UDMI(c,t,0);
	GradTp.x2 = C_UDMI(c,t,1);
	GradTp.x3 = C_UDMI(c,t,2);

	struct dyad2rank d1;
	struct dyad2rank d;
	d1 = dyadicProd(GradTp, GradTp);
	d = dyadMultByConst(d1, con);

	/* the first row of the diadV */
	C_UDSI(c,t,GGx11) = d.x11;
	C_UDSI(c,t,GGx12) = d.x12;
	C_UDSI(c,t,GGx13) = d.x13;

	/* the second row of the diadV */
	C_UDSI(c,t,GGx21) = d.x21;
	C_UDSI(c,t,GGx22) = d.x22;
	C_UDSI(c,t,GGx23) = d.x23;
	
	/* the third row of the diadV */
	C_UDSI(c,t,GGx31) = d.x31;
	C_UDSI(c,t,GGx32) = d.x32;
	C_UDSI(c,t,GGx33) = d.x33;

	struct vector3d divV;
	divV.x1 = C_UDSI_G(c,t,GGx11)[0] + C_UDSI_G(c,t,GGx21)[1] + C_UDSI_G(c,t,GGx31)[2];
	divV.x2 = C_UDSI_G(c,t,GGx12)[0] + C_UDSI_G(c,t,GGx22)[1] + C_UDSI_G(c,t,GGx32)[2];
	divV.x3 = C_UDSI_G(c,t,GGx13)[0] + C_UDSI_G(c,t,GGx23)[1] + C_UDSI_G(c,t,GGx33)[2];

	/* negative due to sign in Eq. (10) */
	C_UDMI(c,t,12) = -divV.x1;
		/*	Message("C_UDMI(c,t,12) = %f%s", C_UDMI(c,t,12), "\n"); */
	C_UDMI(c,t,13) = -divV.x2;
	C_UDMI(c,t,14) = -divV.x3;
}

/* calculates diffusion source terms in the momentum equation */
void diffMomSrcTerms(cell_t c, Thread *t)
{
	double rhos = C_UDMI(c,t,10);
	double rho = C_UDMI(c,t,5);
	double eta = C_UDMI(c,t,7);
	double con = -eta*rhos*C_UDMI(c,t,11)/rho;

	/* values of C_UDSI(c,t,gradGxij) are already calculated in gradTp function 
	   which is run at the begining of ADJUST calcSrcTermsCmpts */
	/* laplcaian calculation (S2 source term components) */
	/* laplacian(G) = div(grad(G)) */
	/* divergence of a vector is calculated as sum of derivatives over columns */
	C_UDMI(c,t,15) = con*(C_UDSI_G(c,t,gradGx11)[0] + C_UDSI_G(c,t,gradGx21)[1] + C_UDSI_G(c,t,gradGx31)[2]);
		/*	Message("C_UDMI(c,t,15) = %f%s", C_UDMI(c,t,15), "\n"); */
	C_UDMI(c,t,16) = con*(C_UDSI_G(c,t,gradGx12)[0] + C_UDSI_G(c,t,gradGx22)[1] + C_UDSI_G(c,t,gradGx32)[2]);
	C_UDMI(c,t,17) = con*(C_UDSI_G(c,t,gradGx13)[0] + C_UDSI_G(c,t,gradGx23)[1] + C_UDSI_G(c,t,gradGx33)[2]);

	/* 1/3-term calculation (S3 source term components) */
	C_UDSI(c,t,divG) = C_UDSI_G(c,t,Gx1)[0] + C_UDSI_G(c,t,Gx2)[1] + C_UDSI_G(c,t,Gx3)[2];
	C_UDMI(c,t,18) = 1./3*con*C_UDSI_G(c,t,divG)[0];
	C_UDMI(c,t,19) = 1./3*con*C_UDSI_G(c,t,divG)[1];
	C_UDMI(c,t,20) = 1./3*con*C_UDSI_G(c,t,divG)[2];
}

/* allows variable storage allocation for UDS */
void uds_derivatives(Domain *d, int n)
{
   /* Code to compute derivative of a variable.  Variable storage allocation first.... */
        MD_Alloc_Storage_Vars(d, SV_UDSI_RG(n), SV_UDSI_G(n), SV_NULL);
        Scalar_Reconstruction(d, SV_UDS_I(n), -1, SV_UDSI_RG(n), NULL);
        Scalar_Derivatives(d, SV_UDS_I(n), -1, SV_UDSI_G(n), SV_UDSI_RG(n), NULL);
        return;
}

/* calculates components of the source terms */
DEFINE_ADJUST(calcSrcTermsCmpts, domain)
{
	/* Make sure there are enough user-defined scalars. */
	if (n_uds < N_REQUIRED_UDS)
		Internal_Error("not enough user-defined scalars allocated");
	
	struct vector3d GradTp;
	double magGradTp;
	Thread *t;
	cell_t c;
   		     
	/* memory allocation for the user defined gradients */
/*   int n;
   for(n=0; n<n_uds; ++n) uds_derivatives(domain, n); */

   /* The UDS gradients are now available */

   {
      Alloc_Storage_Vars(domain, SV_T_RG, SV_T_G,  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_T_RG, SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_P_RG, SV_P_G,  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_P_RG, SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(Gx1), SV_UDSI_G(Gx1),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(Gx1), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(Gx2), SV_UDSI_G(Gx2),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(Gx2), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(Gx3), SV_UDSI_G(Gx3),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(Gx3), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(divG), SV_UDSI_G(divG),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(divG), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(GGx11), SV_UDSI_G(GGx11),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(GGx11), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(GGx12), SV_UDSI_G(GGx12),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(GGx12), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(GGx13), SV_UDSI_G(GGx13),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(GGx13), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(GGx21), SV_UDSI_G(GGx21),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(GGx21), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(GGx22), SV_UDSI_G(GGx22),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(GGx22), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(GGx23), SV_UDSI_G(GGx23),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(GGx23), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(GGx31), SV_UDSI_G(GGx31),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(GGx31), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(GGx32), SV_UDSI_G(GGx32),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(GGx32), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(GGx33), SV_UDSI_G(GGx33),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(GGx33), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(gradGx11), SV_UDSI_G(gradGx11),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(gradGx11), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(gradGx12), SV_UDSI_G(gradGx12),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(gradGx12), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(gradGx13), SV_UDSI_G(gradGx13),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(gradGx13), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(gradGx21), SV_UDSI_G(gradGx21),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(gradGx21), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(gradGx22), SV_UDSI_G(gradGx22),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(gradGx22), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(gradGx23), SV_UDSI_G(gradGx23),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(gradGx23), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(gradGx31), SV_UDSI_G(gradGx31),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(gradGx31), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(gradGx32), SV_UDSI_G(gradGx32),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(gradGx32), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(gradGx33), SV_UDSI_G(gradGx33),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(gradGx33), SV_NULL);
    }



   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(v0e), SV_UDSI_G(v0e),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(v0e), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(v1e), SV_UDSI_G(v1e),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(v1e), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(v2e), SV_UDSI_G(v2e),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(v2e), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(vx1x2), SV_UDSI_G(vx1x2),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(vx1x2), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(vx1y2), SV_UDSI_G(vx1y2),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(vx1y2), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(vx1z2), SV_UDSI_G(vx1z2),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(vx1z2), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(vy1x2), SV_UDSI_G(vy1x2),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(vy1x2), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(vy1y2), SV_UDSI_G(vy1y2),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(vy1y2), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(vy1z2), SV_UDSI_G(vy1z2),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(vy1z2), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(vz1x2), SV_UDSI_G(vz1x2),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(vz1x2), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(vz1y2), SV_UDSI_G(vz1y2),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(vz1y2), SV_NULL);
    }
   {
      Alloc_Storage_Vars(domain, SV_UDSI_RG(vz1z2), SV_UDSI_G(vz1z2),  SV_NULL);
      T_derivatives(domain);
      Free_Storage_Vars(domain, SV_UDSI_RG(vz1z2), SV_NULL);
    }
  
	/* beging loop over all cells */
    thread_loop_c (t,domain)
    {
	    if (
			NULL != THREAD_STORAGE(t,SV_UDS_I(Gx1)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(Gx2)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(Gx3)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(GGx11)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(GGx12)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(GGx13)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(GGx21)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(GGx22)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(GGx23)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(GGx31)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(GGx32)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(GGx33)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(gradGx11)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(gradGx12)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(gradGx13)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(gradGx21)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(gradGx22)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(gradGx23)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(gradGx31)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(gradGx32)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(gradGx33)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(divG)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(v0e)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(v1e)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(v2e)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(vx1x2)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(vx1y2)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(vx1z2)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(vy1x2)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(vy1y2)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(vy1z2)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(vz1x2)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(vz1y2)) &&
		    NULL != THREAD_STORAGE(t,SV_UDS_I(vz1z2)) 
			)  
			{  
			Message("Begin cells loop...\n");
				begin_c_loop (c,t)
				{
					/* Set all UDMs to zero */
					/*for (int i = 0; i < n_udm; i++)
					{
						C_UDMI(c, t, i) = 0.0;
					}*/

					/* calculates G vector */
					GradTp = gradTp(c,t);
					magGradTp = magVector3d(GradTp);

					if (magGradTp < MAG_GRAD_TP_CONV_CRIT) magGradTp = MAG_GRAD_TP_VALUE; 

					/* stores G vector components and |G| */
					C_UDMI(c,t,0) = GradTp.x1; 
					C_UDMI(c,t,1) = GradTp.x2; 
					C_UDMI(c,t,2) = GradTp.x3; 
					C_UDMI(c,t,3) = magGradTp; 

					/* updates k */
					thermCond(c,t); 

				/*	C_UDMI(c,t,5) = rhoArray(c,t); */
					C_UDMI(c,t,8) = AArray(c,t);
					/* calculates and stores rhon in C_UDMI(c,t,9) */
					rhoN(c,t);
					/* calculates and stores rhos in C_UDMI(c,t,10) */
					C_UDMI(c,t,10) = C_UDMI(c,t,5) - C_UDMI(c,t,9);
			
					/* (s/A/rho_n/magGradTp^2)^(1/3) = B */
					C_UDMI(c,t,11) = pow(C_UDMI(c,t,6)/C_UDMI(c,t,8)/C_UDMI(c,t,9)/C_UDMI(c,t,3)/C_UDMI(c,t,3), 1./3);

					/* rho_s/rho*(s/A/rho_n/magGradTp^2)^(1/3)*gradTp_z */
					/* it is used for velocity BC */
					C_UDMI(c,t,21) = C_UDMI(c,t,10)/C_UDMI(c,t,5)*C_UDMI(c,t,11)*C_UDMI(c,t,2);

					/* DOTAD WYDAJE SIE OK */
					/* calculates divergence source term in the momentum equation */
					divDyadV(c,t);

					/* calculates diffusion source terms in the momentum equation */
					diffMomSrcTerms(c,t);     
				}
	        	end_c_loop (c,t)
			}  
    }
/*	Free_Storage_Vars(domain, SV_T_G, SV_NULL);
	Free_Storage_Vars(domain, SV_P_G, SV_NULL); */
}

real rhofN(face_t f, Thread *t, real rho)
{
	const double T = F_T(f,t);
	double rhoNorm = rho*pow(T/TLAMBDA, 5.6);
/*	C_UDMI(c,t,9) = rhon;  */
	return rhoNorm;
}

void rhoN(cell_t c, Thread *t)
{
	double T = C_T(c,t);
	double rhon = C_UDMI(c,t,5)*pow(T/TLAMBDA, 5.6);
	C_UDMI(c,t,9) = rhon;
}

/* calculates thermal conductivity and saves it into CUDMI(c,t,4) */
void thermCond(cell_t c, Thread *t)
{
	/* assigns to magGradTP the magnitude of (-gradP/rho/s + gradT) */
	double magGradTP, k;
	magGradTP = C_UDMI(c,t,3); /* already calculated in calcSrcTermsCmpts in ADJUST */
	if (magGradTP < MAG_GRAD_TP_CONV_CRIT) magGradTP = MAG_GRAD_TP_VALUE;

	k = pow(onebyfArray(c,t)/magGradTP/magGradTP, 1./3.);
	/* stores k */
	C_UDMI(c,t,4) = k;
}

DEFINE_PROPERTY(helThermCond, c, t)
{
	thermCond(c,t);
	double k = C_UDMI(c,t,4);
	return k;
}

DEFINE_PROPERTY(helThermalExp, c, t)
{
	double beta = betaArray(c,t);
	return beta;
}

DEFINE_PROPERTY(helDensity, c, t)
{
	double rho = rhoArray(c,t);
	C_UDMI(c,t,5) = rho;
	return rho;
}

DEFINE_PROPERTY(helDynViscosity, c, t)
{
	double eta = etaArray(c,t);
	C_UDMI(c,t,7) = eta;
	return eta;
}

/*tutaj trzeba uzyc DEFINE_SPECIFIC_HEAT*/
DEFINE_SPECIFIC_HEAT(helSpecificHeat, T, TRef, h, y)
{
	double cp = cpArray(T);
	TRef = TREF;
	*h = cp*(T - TRef);
	return cp;
}

DEFINE_ON_DEMAND(ilehelCond)
{
	Domain *domain;
	domain = Get_Domain(1);
	cell_t c;
	Thread *c_thread;  
	/* c_thread = Lookup_Thread(domain, 12); */
	thread_loop_c(c_thread, domain) /*loops over all cell threads in domain*/
	{
		begin_c_loop(c, c_thread)    /* loops over cells in a cell thread  */
		{
			Message("k = %f%s", C_K_L(c, c_thread), "\n");
		}                         
		end_c_loop(c, c_thread)
	}
}

DEFINE_ON_DEMAND(sprDefine)
{
	Domain *domain;
	domain = Get_Domain(1);
	cell_t c;
	Thread *c_thread;  
	/* c_thread = Lookup_Thread(domain, 12); */
	thread_loop_c(c_thread, domain) /*loops over all cell threads in domain*/
	{
		begin_c_loop(c, c_thread)    /* loops over cells in a cell thread  */
		{
			Message("Tlambda = %f%s", TLAMBDA, "\n");
			Message("MAG_GRAD_TP_CONV_CRIT = %f%s", MAG_GRAD_TP_CONV_CRIT, "\n");
			Message("MAG_GRAD_TP_VALUE = %f%s", MAG_GRAD_TP_CONV_CRIT, "\n");
			Message("n_uds = %d%s", n_uds, "\n"); 
		}                         
		end_c_loop(c, c_thread)
	}
}

DEFINE_ADJUST(sprTabWlTermHelu, domain)
{
	cell_t c;
	Thread *t;  
	/* c_thread = Lookup_Thread(domain, 12); */
	thread_loop_c(t, domain) /*loops over all cell threads in domain*/
	{
		begin_c_loop(c, t)    /* loops over cells in a cell thread  */
		{
			double T = C_T(c,t);
			Message("rho = %f%s", rhoArray(c,t), "\n");
			Message("cp = %f%s", cpArray(T), "\n");
			Message("eta = %f%s", etaArray(c,t), "\n");
			Message("s = %f%s", sArray(c,t), "\n");
			Message("A = %f%s", AArray(c,t), "\n");
			Message("beta = %f%s", betaArray(c,t), "\n");
			Message("1/f = %f%s", onebyfArray(c,t), "\n"); 
		}                         
		end_c_loop(c, t)
	}
}

/* internal convection source term in the energy equation */
DEFINE_SOURCE(energyInternalConvSrcTerm, c, t, dS, eqn)
{
	/* internal convection term */
	double k = C_UDMI(c,t,4);
	double dvdx, dvdy, dvdz, divkGradTp, rho, s;

	rho = C_UDMI(c,t,5);
	s = C_UDMI(c,t,6);

	C_UDSI(c,t,v0e) = -k*C_P_G(c,t)[0]/rho/s;
	C_UDSI(c,t,v1e) = -k*C_P_G(c,t)[1]/rho/s;
	C_UDSI(c,t,v2e) = -k*C_P_G(c,t)[2]/rho/s;

	dvdx = C_UDSI_G(c,t,v0e)[0];
	dvdy = C_UDSI_G(c,t,v1e)[1];
	dvdz = C_UDSI_G(c,t,v2e)[2];

	/* jak bedzie ciezko ze zbieznoscia to trzeba zrobic implicit */
	divkGradTp = dvdx + dvdy + dvdz;

	dS[eqn] = 0.0;

	return divkGradTp;
}


/* Gorter-Mellink source term in the energy equation */
DEFINE_SOURCE(energyGorterMellinkiSrcTerm, c, t, dS, eqn)
{
	/* Gorter-Mellink mutual friction term (GM) */
	double GM, A, s, rhon, rhos, magGradTp;
	A = C_UDMI(c,t,8);
	s = C_UDMI(c,t,6);
	rhon = C_UDMI(c,t,9);
	rhos = C_UDMI(c,t,10);
	magGradTp = C_UDMI(c,t,3);

	GM = A*rhon*rhos*pow(C_UDMI(c,t,11)*C_UDMI(c,t,3), 4.0);

	/* jak bedzie ciezko ze zbieznoscia to trzeba zrobic implicit */
	dS[eqn] = 0.0;

	return GM;
}

/* divergence source terms in the x momentum equation */
DEFINE_SOURCE(divMom_x_SrcTerm, c, t, dS, eqn)
{
	double divmomxsrcterm = C_UDMI(c,t,12); 

	dS[eqn] = 0.0;
	return divmomxsrcterm;
}

/* divergence source terms in the y momentum equation */
DEFINE_SOURCE(divMom_y_SrcTerm, c, t, dS, eqn)
{
	double divmomysrcterm = C_UDMI(c,t,13); 

	dS[eqn] = 0.0;
	return divmomysrcterm;
}

/* divergence source term in the z momentum equation */
DEFINE_SOURCE(divMom_z_SrcTerm, c, t, dS, eqn)
{
	double divmomzsrcterm = C_UDMI(c,t,14); 

	dS[eqn] = 0.0;
	return divmomzsrcterm;
}

/* laplacian source term in the x momentum equation */
DEFINE_SOURCE(laplMom_x_SrcTerm, c, t, dS, eqn)
{
	double laplmomxsrcterm = C_UDMI(c,t,15); 

	dS[eqn] = 0.0;
	return laplmomxsrcterm;
}

/* laplacian source term in the y momentum equation */
DEFINE_SOURCE(laplMom_y_SrcTerm, c, t, dS, eqn)
{
	double laplmomysrcterm = C_UDMI(c,t,16); 

	dS[eqn] = 0.0;
	return laplmomysrcterm;
}

/* laplacian source term in the z momentum equation */
DEFINE_SOURCE(laplMom_z_SrcTerm, c, t, dS, eqn)
{
	double laplmomzsrcterm = C_UDMI(c,t,17); 

	dS[eqn] = 0.0;
	return laplmomzsrcterm;
}

/* 1/3-term source term in the x momentum equation */
DEFINE_SOURCE(oneThirdMom_x_SrcTerm, c, t, dS, eqn)
{
	double onethirdmomxsrcterm = C_UDMI(c,t,18); 

	dS[eqn] = 0.0;
	return onethirdmomxsrcterm;
}

/* 1/3-term source term in the y momentum equation */
DEFINE_SOURCE(oneThirdMom_y_SrcTerm, c, t, dS, eqn)
{
	double onethirdmomysrcterm = C_UDMI(c,t,19); 

	dS[eqn] = 0.0;
	return onethirdmomysrcterm;
}
/* 1/3-term source term in the z momentum equation */
DEFINE_SOURCE(oneThirdMom_z_SrcTerm, c, t, dS, eqn)
{
	double onethirdmomzsrcterm = C_UDMI(c,t,20); 

	dS[eqn] = 0.0;
	return onethirdmomzsrcterm;
}

DEFINE_PROFILE(wallVelocity_z,tf,i)
{
	cell_t c0;
	face_t f;
	Thread *t0;

	begin_f_loop(f,tf)
		c0 = F_C0(f,tf);
		t0 = THREAD_T0(tf);
		F_PROFILE(f,tf,i) = C_UDMI(c0,t0,21);
	end_f_loop(f,tf)
}

/* initial temperature at the begiging of calculations in domain - SP */

DEFINE_INIT(my_init_temp,domain)
{
	cell_t c;
	Thread *ct; 
	double xc[ND_ND], z;
	const double Tb = 1.80; /* temperature at the inlet */
	const double Te = 1.70; /* expected temperature at the outlet */
	const double hd = 0.0058; /* lenght of the domain in z */
	
	thread_loop_c(ct, domain)
	{
		begin_c_loop(c,ct)
		{
			C_CENTROID(xc,c,ct);
			z = xc[2];
			C_T(c,ct) = Tb - (((Te-Tb)/hd)*z);
		}
		end_c_loop(c,ct)
	}
}
