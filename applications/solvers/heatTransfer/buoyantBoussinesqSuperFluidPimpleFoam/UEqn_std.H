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
        rhon*rhos/rho*B*B*G*G
//        mesh,
//        dimensionedTensor("V", dimensionSet(0,2,-2,0,0,0,0), tensor(0))
    );
//    V = G*G*rhon*rhos*B*B/rho/rho;

    // source term M
    volVectorField M
    (
        IOobject
        (
            "M",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(V)/rho
        + turbulence->nuEff()*rhos/rho*B*
    //  (4./3*fvc::laplacian(G) + 1./3*fvc::curl(fvc::curl(G)))
          (fvc::laplacian(G) + 1./3*fvc::grad(fvc::div(G)))
    );
    Info<< "TUTAJ" << endl;

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
      - M
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
    //          - (fvc::interpolate(M) & mesh.Sf())
            )
        );

        fvOptions.correct(U);
    }
