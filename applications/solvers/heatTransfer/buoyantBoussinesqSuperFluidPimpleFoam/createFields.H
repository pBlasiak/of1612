Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Dla test case wycinek
//const dimensionedScalar TMin("TMin", dimensionSet(0,0,0,1,0,0,0), 1.675);
//const dimensionedScalar TMax("TMax", dimensionSet(0,0,0,1,0,0,0), 1.7);
//T = (TMax - TMin)/(max(mesh.C().component(0)) - min(mesh.C().component(0)))*
//        mesh.C().component(0) + TMin;

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//Info<< "Reading field p_rgh\n" << endl;
//volScalarField p_rgh
//(
//    IOobject
//    (
//        "p_rgh",
//        runTime.timeName(),
//        mesh,
//        IOobject::MUST_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh
//);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Creating simplifiedSuperFluid model\n" << endl;
autoPtr<simplifiedSuperFluid> superFluid =
    simplifiedSuperFluid::New(U, phi);

superFluid->correct();

#include "readTransportProperties.H"

volScalarField cp
(
    IOobject
    (
        "cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	superFluid->cpHe()
);

Info<< "Creating turbulence model\n" << endl;
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, superFluid())
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	superFluid->rhoHe()
);

volScalarField Pr
(
    IOobject
    (
        "Pr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	superFluid->Pr()
);

// Thermal expansion coefficient [1/K]
volScalarField beta
(
    IOobject
    (
        "beta",
        runTime.timeName(),
        mesh
    ),
    superFluid->betaHe()
);

// Kinematic density for buoyancy force
volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh
    ),
    1.0 - beta*(T - TRef)
);

// kinematic turbulent thermal thermal conductivity m2/s
Info<< "Reading field alphat\n" << endl;
volScalarField alphat
(
    IOobject
    (
        "alphat",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField alphaEff
(
    IOobject
    (
        "alphaEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    superFluid->alpha()/rho/cp
);

//volScalarField nut
//(
//    IOobject
//    (
//        "nut",
//        runTime.timeName(),
//        mesh,
//        IOobject::MUST_READ,
//        IOobject::AUTO_WRITE
//    ),
//	mesh
////	dimensionedScalar("nut", dimensionSet(0,2,-1,0,0,0,0), scalar(0))
////    turbulence->nut()()
//);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


//volScalarField p
//(
//    IOobject
//    (
//        "p",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::AUTO_WRITE
//    ),
//    p_rgh + rhok*gh
//);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
//    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

//if (p_rgh.needReference())
//{
//    p += dimensionedScalar
//    (
//        "p",
//        p.dimensions(),
//        pRefValue - getRefCellValue(p, pRefCell)
//    );
//}

//mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(p.name());

// normal velocity component
volVectorField Un
(
    IOobject
    (
        "Un",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// super velocity component
volVectorField Us
(
    IOobject
    (
        "Us",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField UnmUs
(
    IOobject
    (
        "UnmUs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Un - Us
);

// G = -gradp/rho/s + gradT
volVectorField G
(
    IOobject
    (
        "G",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	superFluid->G()
//    mesh,
//    dimensionedVector("G", dimensionSet(0,-1,0,1,0,0,0), vector(0,0,0))
);

volScalarField rhon
(
    IOobject
    (
        "rhon",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	superFluid->rhon()
);
Info<< "rhon = " << rhon << endl;

volScalarField rhos
(
    IOobject
    (
        "rhos",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	superFluid->rhos()
);
Info<< "rhos = " << rhos << endl;

volScalarField s
(
    IOobject
    (
        "s",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	superFluid->sHe()
);

// B = (s/A/rhon/magG^2)^(1./3)
volScalarField B
(
    IOobject
    (
        "B",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	superFluid->B()
//    dimensionedScalar("B", dimensionSet(0,2,-1,-1,0,0,0), scalar(0))
);

surfaceScalarField q
(
    IOobject
    (
        "q",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(superFluid->alpha())*fvc::snGrad(T)
);

