    // Solve the momentum equation

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

//    fvVectorMatrix UEqn
//    (
//        fvm::ddt(U) + fvm::div(phi, U)
//      + MRF.DDt(U)
//      + turbulence->divDevReff(U)
//     ==
//        fvOptions(U)
//    );

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
            )
        );

        fvOptions.correct(U);
    }
