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

#include "Helium.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(Helium, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        Helium,
        dictionary
    );

	const Foam::dimensionedScalar
	Foam::viscosityModels::Helium::Tlambda_("Tlambda", dimTemperature, 2.1711132461);

	const Foam::dimensionedScalar
	Foam::viscosityModels::Helium::TMin_("TMin", dimTemperature, 1.5);
	
	const Foam::dimensionedScalar
	Foam::viscosityModels::Helium::TMax_("TMax", dimTemperature, 2.167);
	
	const Foam::label
	Foam::viscosityModels::Helium::indexMin_(0);
	
	const Foam::label
	Foam::viscosityModels::Helium::indexMax_(667);
	
	const Foam::scalar
	Foam::viscosityModels::Helium::dT_(0.001);

	#include "staticTables.H"
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::Helium::calcNu() 
{
	calcHeProp(etaHe_, etaHeTable_);

    volScalarField nu
    (
        IOobject
        (
            "nuLocal",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
		etaHe_/rhoHe_
    );

    return max
    (
        nuMin_,
        min
        (
            nuMax_,
			nu
        )
    );
}

void Foam::viscosityModels::Helium::calcHeProp
(
	Foam::volScalarField& vsf,
	const List<scalar>& vsfTable
)
{
	forAll(T_, celli)
	{
		if (T_[celli] < TMin_.value())
		{
			vsf[celli] = vsfTable[indexMin_];
		}
		else if (T_[celli] > TMax_.value())
		{
			vsf[celli] = vsfTable[indexMax_];
		}
		else
		{
			label index = (T_[celli] - TMin_.value())/dT_;
			scalar Ti1 = TMin_.value() + index*dT_;
			scalar Ti2 = Ti1 + dT_;
			scalar a = (vsfTable[index + 1] - vsfTable[index])/(Ti2 - Ti1);
			scalar b = vsfTable[index] - a*Ti1;
			vsf[celli] =  a*T_[celli] + b;
		}
	}

	forAll(T_.boundaryField(), patchi)
	{
		forAll(T_.boundaryField()[patchi], i)
		{
			if (T_.boundaryField()[patchi][i] < TMin_.value())
			{
				vsf.boundaryFieldRef()[patchi][i] = vsfTable[indexMin_];
			}
			else if (T_.boundaryField()[patchi][i] > TMax_.value())
			{
				vsf.boundaryFieldRef()[patchi][i] = vsfTable[indexMax_];
			}
			else
			{
				label index = (T_.boundaryField()[patchi][i] - TMin_.value())/dT_;
				scalar Ti1 = TMin_.value() + index*dT_;
				scalar Ti2 = Ti1 + dT_;
				scalar a = (vsfTable[index + 1] - vsfTable[index])/(Ti2 - Ti1);
				scalar b = vsfTable[index] - a*Ti1;
				vsf.boundaryFieldRef()[patchi][i] =  a*T_.boundaryField()[patchi][i] + b;
			}
		}
	}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::Helium::Helium
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    HeliumCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
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
		dimensionedScalar("betaHe", dimless/dimTemperature, 0.0)
    ),

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
		dimensionedScalar("AGM", dimensionSet(-1,1,1,0,0,0,0), 0.0)
    ),

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
		dimensionedScalar("sHe", dimensionSet(0,2,-2,-1,0,0,0), 0.0)
    ),

    etaHe_
    (
        IOobject
        (
            "etaHe",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
		dimensionedScalar("etaHe", dimensionSet(1,-1,-1,0,0,0,0), 0.0)
    ),

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
		dimensionedScalar("cpHe", dimensionSet(0,2,-2,-1,0,0,0), 0.0)
    ),
	

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
		dimensionedScalar("onebyf", dimensionSet(3,1,-9,-1,0,0,0), 0.0)
    ),
    rhoHe_("rhoHe", dimDensity, HeliumCoeffs_),
	T_(U_.db().lookupObject<volScalarField>("T")),
    nuMin_("nuMin", dimViscosity, etaHeTable_[indexMin_]/rhoHe_.value()),
    nuMax_("nuMax", dimViscosity, etaHeTable_[indexMax_]/rhoHe_.value()),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{
	Info<< "rhoHe = " << rhoHe_ << endl;
	Info<< "nuMin = " << nuMin_ << endl;
	Info<< "nuMax = " << nuMax_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::viscosityModels::Helium::correct()
{
	Info<< "Helium updates thermal properties..." << endl;
	nu_ = calcNu();
	calcHeProp(betaHe_, betaHeTable_);
	calcHeProp(AGMHe_, AGMHeTable_);
	calcHeProp(sHe_, sHeTable_);
	calcHeProp(cpHe_, cpHeTable_);
	calcHeProp(onebyf_, onebyfTable_);
}

bool Foam::viscosityModels::Helium::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    HeliumCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    HeliumCoeffs_.lookup("rhoHe") >> rhoHe_;

    return true;
}


// ************************************************************************* //
