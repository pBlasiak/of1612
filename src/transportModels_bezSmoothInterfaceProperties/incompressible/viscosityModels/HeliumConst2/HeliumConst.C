/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "HeliumConst.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(HeliumConst, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        HeliumConst,
        dictionary
    );

	const Foam::dimensionedScalar
	Foam::viscosityModels::HeliumConst::Tlambda_("Tlambda", dimTemperature, 2.1711132461);
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::HeliumConst::HeliumConst
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    HeliumConstCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    betaHe0_("betaHe", dimless/dimTemperature, HeliumConstCoeffs_),
    betaHe_
    (
        IOobject
        (
            "betaHe",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
		betaHe0_
    ),

    AGMHe0_("AGMHe", dimensionSet(-1,1,1,0,0,0,0), HeliumConstCoeffs_),
    AGMHe_
    (
        IOobject
        (
            "AGMHe",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
		AGMHe0_
    ),

    sHe0_("sHe", dimensionSet(0,2,-2,-1,0,0,0), HeliumConstCoeffs_),
    sHe_
    (
        IOobject
        (
            "sHe",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
		sHe0_
    ),

    cpHe0_("cpHe", dimensionSet(0,2,-2,-1,0,0,0), HeliumConstCoeffs_),
    cpHe_
    (
        IOobject
        (
            "cpHe",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
		cpHe0_
    ),
	
    onebyf0_("onebyf", dimensionSet(3,1,-9,-1,0,0,0), HeliumConstCoeffs_),
    onebyf_
    (
        IOobject
        (
            "onebyf",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
		onebyf0_
    ),
    rhoHe_("rhoHe", dimDensity, HeliumConstCoeffs_),
    nu0_("nu", dimViscosity, HeliumConstCoeffs_),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh(),
        nu0_
    )
{
	Info<< "betaHe = " << betaHe_ << endl;
	Info<< "AGMHe = " << AGMHe_ << endl;
	Info<< "sHe = " << sHe_ << endl;
	Info<< "cpHe = " << cpHe_ << endl;
	Info<< "onebyf = " << onebyf_ << endl;
	Info<< "rhoHe = " << rhoHe_ << endl;
	Info<< "nuHe = " << nu_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::HeliumConst::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    HeliumConstCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    HeliumConstCoeffs_.lookup("betaHe") >> betaHe0_;
	betaHe_ = betaHe0_;
    HeliumConstCoeffs_.lookup("AGMHe") >> AGMHe0_;
	AGMHe_ = AGMHe0_;
    HeliumConstCoeffs_.lookup("sHe") >> sHe0_;
	sHe_ = sHe0_;
    HeliumConstCoeffs_.lookup("cpHe") >> cpHe0_;
	cpHe_ = cpHe0_;
    HeliumConstCoeffs_.lookup("onebyf") >> onebyf0_;
	onebyf_ = onebyf0_;
    HeliumConstCoeffs_.lookup("rhoHe") >> rhoHe_;
    HeliumConstCoeffs_.lookup("nu") >> nu0_;
    nu_ = nu0_;
	Info<< "betaHe = " << betaHe_ << endl;
	Info<< "AGMHe = "  << AGMHe_  << endl;
	Info<< "sHe = "    << sHe_    << endl;
	Info<< "cpHe = "   << cpHe_   << endl;
	Info<< "onebyf = " << onebyf_ << endl;
	Info<< "rhoHe = "  << rhoHe_  << endl;
	Info<< "nuHe = "   << nu_     << endl;

    return true;
}


// ************************************************************************* //
