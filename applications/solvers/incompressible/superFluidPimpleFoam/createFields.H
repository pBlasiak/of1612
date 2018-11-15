Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

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

Info<< "Reading field T\n" << endl;
const dimensionedScalar TMin("TMin", dimensionSet(0,0,0,1,0,0,0), 1.7);
const dimensionedScalar TMax("TMax", dimensionSet(0,0,0,1,0,0,0), 1.8);
                         ;
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
T = (TMax - TMin)/(max(mesh.C().component(2)) - min(mesh.C().component(2)))*
        mesh.C().component(2) + TMin;


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


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

#include "readTransportProperties.H"

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

const dimensionedScalar AGM("AGM", dimensionSet(-1,1,1,0,0,0,0), 764.3353);
const dimensionedScalar s("s", dimensionSet(0,2,-2,-1,0,0,0), 401.43055);
const dimensionedScalar cp("cp", dimensionSet(0,2,-2,-1,0,0,0), 2207.427);
const dimensionedScalar rho("rho", dimensionSet(1,-3,0,0,0,0,0), 147.0844);
const dimensionedScalar onebyf("onebyf", dimensionSet(3,1,-9,-1,0,0,0), 6445631951700.62);
const dimensionedScalar Tlambda("Tlambda", dimensionSet(0,0,0,1,0,0,0), 2.1768);

dimensionedScalar MAG_G_CONV_CRIT
(
    "MAG_G_CONV_CRIT",
    dimless,
    transportProperties.lookup("MAG_G_CONV_CRIT")
);

dimensionedScalar MAG_G_MIN_VALUE
(
    "MAG_G_MIN_VALUE",
    dimless,
    transportProperties.lookup("MAG_G_MIN_VALUE")
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
    mesh,
    dimensionedVector("G", dimensionSet(0,-1,0,1,0,0,0), vector(0,0,0))
);

volScalarField magG
(
    IOobject
    (
        "magG",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("magG", dimensionSet(0,-1,0,1,0,0,0), 1e-3)
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
    mesh,
    rho
);

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
    mesh,
    rho
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
    mesh,
    dimensionedScalar("B", dimensionSet(0,2,-1,-1,0,0,0), scalar(0))
);
