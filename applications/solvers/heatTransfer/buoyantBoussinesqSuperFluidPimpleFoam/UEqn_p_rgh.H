    // Solve the momentum equation

    MRF.correctBoundaryVelocity(U);

//    // diad V
//    volTensorField V
//    (
//        IOobject
//        (
//            "V",
//            runTime.timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        mesh,
//        dimensionedTensor("V", dimensionSet(0,2,-2,0,0,0,0), tensor(0))
//    );
//    V = G*G*rhon*rhos*B*B/rho/rho;


	nut = turbulence->nut();

    // source term M
    tmp<volVectorField> tM
    (
		superFluid->M1()/rho
      + (superFluid->nu() + nut)*superFluid->M2()
    );
	volVectorField& M = tM.ref();

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
              - (fvc::interpolate(M) & mesh.Sf())
            )
        );

        fvOptions.correct(U);
    }
