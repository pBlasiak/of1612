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

	const Foam::dimensionedScalar
	Foam::viscosityModels::HeliumConst::TMin_("TMin", dimTemperature, 1.5);
	
	const Foam::dimensionedScalar
	Foam::viscosityModels::HeliumConst::TMax_("TMax", dimTemperature, 2.167);
	
	const Foam::label
	Foam::viscosityModels::HeliumConst::indexMin_(0);
	
	const Foam::label
	Foam::viscosityModels::HeliumConst::indexMax_(667);
	
	const Foam::scalar
	Foam::viscosityModels::HeliumConst::dT_(0.001);

	#include "staticTables.H"
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::HeliumConst::calcNu() 
{
	calcHeProp(etaHe_, etaHeTable_);

    tmp<volScalarField> nu
    (
		new volScalarField
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
		)
    );

    return nu;
}

void Foam::viscosityModels::HeliumConst::calcHeProp
(
	Foam::volScalarField& vsf,
	const List<scalar>& vsfTable
)
{
	forAll(vsf, celli)
	{
		if (TMean_ < TMin_)
		{
			vsf[celli] = vsfTable[indexMin_];
		}
		else if (TMean_ > TMax_)
		{
			vsf[celli] = vsfTable[indexMax_];
		}
		else
		{
			label index = (TMean_.value() - TMin_.value())/dT_;
			scalar Ti1 = TMin_.value() + index*dT_;
			scalar Ti2 = Ti1 + dT_;
			scalar a = (vsfTable[index + 1] - vsfTable[index])/(Ti2 - Ti1);
			scalar b = vsfTable[index] - a*Ti1;
			scalar value = a*TMean_.value() + b;
			vsf[celli] = value;
		}
	}

	forAll(vsf.boundaryField(), patchi)
	{
		forAll(vsf.boundaryField()[patchi], i)
		{
			if (TMean_ < TMin_)
			{
				vsf.boundaryFieldRef()[patchi][i] = vsfTable[indexMin_];
			}
			else if (TMean_ > TMax_)
			{
				vsf.boundaryFieldRef()[patchi][i] = vsfTable[indexMax_];
			}
			else
			{
				label index = (TMean_.value() - TMin_.value())/dT_;
				scalar Ti1 = TMin_.value() + index*dT_;
				scalar Ti2 = Ti1 + dT_;
				scalar a = (vsfTable[index + 1] - vsfTable[index])/(Ti2 - Ti1);
				scalar b = vsfTable[index] - a*Ti1;
				scalar value = a*TMean_.value() + b;
				vsf.boundaryFieldRef()[patchi][i] = value;
			}
		}
	}
}

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
    TMean_("TMean", dimTemperature, HeliumConstCoeffs_),
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
    rhoHe_("rhoHe", dimDensity, HeliumConstCoeffs_),
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
	nu_ = calcNu();
	calcHeProp(betaHe_, betaHeTable_);
	calcHeProp(AGMHe_, AGMHeTable_);
	calcHeProp(sHe_, sHeTable_);
	calcHeProp(cpHe_, cpHeTable_);
	calcHeProp(onebyf_, onebyfTable_);

	Info<< "TMean = " << TMean_ << endl;
	Info<< "betaHe = " << betaHe_ << endl;
	Info<< "AGMHe = " << AGMHe_ << endl;
	Info<< "sHe = " << sHe_ << endl;
	Info<< "cpHe = " << cpHe_ << endl;
	Info<< "onebyf = " << onebyf_ << endl;
	Info<< "rhoHe = " << rhoHe_ << endl;
	Info<< "nuHe = " << nu_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::viscosityModels::HeliumConst::correct()
{}

bool Foam::viscosityModels::HeliumConst::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    HeliumConstCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    HeliumConstCoeffs_.lookup("rhoHe") >> rhoHe_;
    HeliumConstCoeffs_.lookup("TMean") >> TMean_;

	Info<< "HeliumConst updates thermal properties..." << endl;
	nu_ = calcNu();
	calcHeProp(betaHe_, betaHeTable_);
	calcHeProp(AGMHe_, AGMHeTable_);
	calcHeProp(sHe_, sHeTable_);
	calcHeProp(cpHe_, cpHeTable_);
	calcHeProp(onebyf_, onebyfTable_);

	Info<< "TMean = " << TMean_ << endl;
	Info<< "betaHe = " << betaHe_ << endl;
	Info<< "AGMHe = " << AGMHe_ << endl;
	Info<< "sHe = " << sHe_ << endl;
	Info<< "cpHe = " << cpHe_ << endl;
	Info<< "onebyf = " << onebyf_ << endl;
	Info<< "rhoHe = " << rhoHe_ << endl;
	Info<< "nuHe = " << nu_ << endl;

    return true;
}


// ************************************************************************* //
