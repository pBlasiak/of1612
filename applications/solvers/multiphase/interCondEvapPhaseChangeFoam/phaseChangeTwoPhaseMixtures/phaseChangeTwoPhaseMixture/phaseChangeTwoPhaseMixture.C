/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "phaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseChangeTwoPhaseMixture, 0);
    defineRunTimeSelectionTable(phaseChangeTwoPhaseMixture, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixture::phaseChangeTwoPhaseMixture
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    thermalIncompressibleTwoPhaseMixture(U, phi),
    phaseChangeTwoPhaseMixtureCoeffs_(subDict(type + "Coeffs")),
    TSat_("TSat", dimTemperature, phaseChangeTwoPhaseMixtureCoeffs_.lookup("TSat")),
    pSat_("pSat", dimPressure, phaseChangeTwoPhaseMixtureCoeffs_.lookup("pSat")),
    hEvap_("hEvap", dimEnergy/dimMass, phaseChangeTwoPhaseMixtureCoeffs_.lookup("hEvap")),
    R_("R", dimGasConstant, phaseChangeTwoPhaseMixtureCoeffs_.lookup("R")),
    TSatLocal_(readBool(phaseChangeTwoPhaseMixtureCoeffs_.lookup("TSatLocal"))),
    cutoff_("cutoff", dimless, phaseChangeTwoPhaseMixtureCoeffs_.lookup("cutoff"))
{
	Info<< "pSat = "    << pSat_ << endl;
	Info<< "TSat = "    << TSat_ << endl;
	Info<< "hEvap = "   << hEvap_ << endl;
	Info<< "R = "       << R_ << endl;
	Info<< "cutoff = "  << cutoff_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::volScalarField
Foam::phaseChangeTwoPhaseMixture::TSatLocal() const
{
	if (TSatLocal_)
	{
		Info <<"TSatlocal" << endl;
	    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
	    volScalarField oneByT = 1.0/TSat_ - R_/hEvap_*log(max(p/pSat_,1E-08));

	    return volScalarField
	    (
	           1.0/oneByT
	    );
	}
	else
	{
		Info <<"TSat" << endl;
		volScalarField one
		(
			IOobject
			(
				"one",
				alpha1_.mesh(),
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			alpha1_.mesh(),
			dimensionedScalar ("one",dimless, 1.0)
	    );

		return volScalarField
		(
			one*TSat_
	    );
	}
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixture::vDotAlphal() const
{
    volScalarField alphalCoeff(1.0/rho1() - alpha1_*(1.0/rho1() - 1.0/rho2()));
    Pair<tmp<volScalarField>> mDotAlphal = this->mDotAlphal();

    return Pair<tmp<volScalarField>>
    (
        alphalCoeff*mDotAlphal[0],
        alphalCoeff*mDotAlphal[1]
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixture::vDotP() const
{
    dimensionedScalar pCoeff(1.0/rho1() - 1.0/rho2());
    Pair<tmp<volScalarField>> mDotP = this->mDotP();

    return Pair<tmp<volScalarField>>(pCoeff*mDotP[0], pCoeff*mDotP[1]);
}


bool Foam::phaseChangeTwoPhaseMixture::read()
{
    if (incompressibleTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("TSat") >> TSat_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("pSat") >> pSat_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("cutoff") >> cutoff_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
