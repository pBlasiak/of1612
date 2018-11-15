{
    volScalarField kEff
    (
        IOobject
        (
            "kEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kEff", dimensionSet(0,2,-1,0,0,0,0), scalar(0))
    );
    kEff =
        pow(max(onebyf/magG/magG, dimensionedScalar("small", dimensionSet(3,3,-9,-3,0,0,0), SMALL)), 1./3)/rho/cp
      + turbulence->nut()/Prt;

    surfaceScalarField kfEff
    (
        IOobject
        (
            "kfEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kfEff", dimensionSet(0,2,-1,0,0,0,0), scalar(0))
    );
    kfEff =
        fvc::interpolate(pow(max(onebyf/magG/magG, dimensionedScalar("small", dimensionSet(3,3,-9,-3,0,0,0), SMALL)), 1./3))/rho/cp
      + fvc::interpolate(turbulence->nut())/Prt;

    // Goerter-Mellink source term
    volScalarField GM(AGM*rhon*rhos*pow(max(B*magG, dimensionedScalar("small", dimensionSet(0,1,-1,0,0,0,0), SMALL)), 4.0));
    volScalarField pDivU = p*fvc::div(phi)/cp;

//    //Energy-temperature loop
//    label nEnergyLoops(readLabel(pimple.dict().lookup("nEnergyLoops")));
//    for (int EEqnCount=0; EEqnCount < nEnergyLoops; EEqnCount++)
//    {
    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
//      - fvm::laplacian(kfEff, T)
      - fvm::laplacian(kEff, T)
//      - fvc::laplacian(kEff, T)
//      - kEff*fvm::laplacian(T)
//      - fvm::Sp(fvc::div(kEff*gradp)/T, T)
     ==
        fvm::Sp(fvc::div(kEff*gradp)/T, T)
//        kEff*fvc::div(gradp)
//      - fvc::laplacian(fvc::interpolate(kEff/rho/s), p*rho)
//      - fvc::laplacian(kEff/rho/s, p*rho)
//      - kEff/rho/s*fvc::laplacian(p*rho)
//      + GM/rho/cp
//      + fvm::Sp(GM/rho/cp/T, T)
      + fvm::SuSp(GM/rho/cp/T, T)
//      - fvm::Sp(pDivU/T, T)
      - fvm::SuSp(pDivU/T, T)
//      - pDivU
    );


    TEqn.relax();
    TEqn.solve();
//    }

    T.correctBoundaryConditions();
    Info<< "min/max(T) = " << min(T).value() << ", "
        << max(T).value() <<endl;
}
