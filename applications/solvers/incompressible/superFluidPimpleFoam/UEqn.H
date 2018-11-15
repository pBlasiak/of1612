// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

// diad V
volTensorField V
(
    IOobject
    (
        "V",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor("V", dimensionSet(0,2,-2,0,0,0,0), tensor(0))
);
V = G*G*rhon*rhos*B*B/rho/rho;


tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
  - fvc::div(V)
  - turbulence->nuEff()*rhos/rho*B*
//  (4./3*fvc::laplacian(G) + 1./3*fvc::curl(fvc::curl(G)))
    (fvc::laplacian(G) + 1./3*fvc::grad(fvc::div(G)))
  + fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

// RHS without -grad(p)
//const surfaceScalarField momSource =
//                fvc::interpolate
//                (
//                  - fvc::div(V)
//                  - turbulence->nuEff()*rhos/rho*B*
//                    (4./3*fvc::laplacian(G) + 1./3*fvc::curl(fvc::curl(G)))
//                      (fvc::laplacian(G) + 1./3*fvc::grad(fvc::div(G)))
//                ) & mesh.Sf();

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));
//    solve
//    (
//        UEqn
//     ==
//        fvc::reconstruct
//        (
//            momSource
//          - fvc::snGrad(p)*mesh.magSf()
//        )
//    );

    fvOptions.correct(U);
}
