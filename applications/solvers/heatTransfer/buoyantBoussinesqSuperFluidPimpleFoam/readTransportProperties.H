// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, superFluid());

// Laminar Prandtl number
//dimensionedScalar Pr("Pr", dimless, superFluid());

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, superFluid());

//Info<< "Creating variables AGM, s, cp, rho, onebyf, Tlambda\n" << endl;
//volScalarField AGM
//(
//    IOobject
//    (
//        "AGM",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::NO_WRITE
//    ),
//    mesh,
//    dimensionedScalar("AGM", dimensionSet(-1,1,1,0,0,0,0), superFluid())
//);
//
//volScalarField s
//(
//    IOobject
//    (
//        "s",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::NO_WRITE
//    ),
//    mesh,
//    dimensionedScalar("s", dimensionSet(0,2,-2,-1,0,0,0), superFluid())
//);


//volScalarField onebyf
//(
//    IOobject
//    (
//        "onebyf",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::NO_WRITE
//    ),
//    mesh,
//    dimensionedScalar("onebyf", dimensionSet(3,1,-9,-1,0,0,0), superFluid())
//);
//
//volScalarField Tlambda
//(
//    IOobject
//    (
//        "Tlambda",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::NO_WRITE
//    ),
//    mesh,
//    dimensionedScalar("Tlambda", dimensionSet(0,0,0,1,0,0,0), superFluid())
//);
