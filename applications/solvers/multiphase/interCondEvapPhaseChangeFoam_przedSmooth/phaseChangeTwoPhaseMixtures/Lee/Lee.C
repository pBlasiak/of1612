/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "Lee.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(Lee, 0);
    addToRunTimeSelectionTable(phaseChangeTwoPhaseMixture, Lee, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixtures::Lee::Lee
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    phaseChangeTwoPhaseMixture(typeName, U, phi),

    Cc_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cc")),
    Cv_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cv")),

    p0_("p0", pSat().dimensions(), 0.0),
    TSat_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("TSat")),
    T0_("T0", dimensionSet(0,0,0,0,0,0,0), 0.0),

    mcCoeff_(Cc_*rho2()),
    mvCoeff_(Cv_*rho1())
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Lee::mDot() const
{
    const volScalarField& T = alpha1_.db().lookupObject<volScalarField>("T");
    volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));

    volScalarField mCon
    (
        IOobject
        (
            "mCon",
            U_.time().timeName(),
            U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("mCon", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0)
    );
    
    volScalarField mVap
    (
        IOobject
        (
            "mVap",
            U_.time().timeName(),
            U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("mVap", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0)
    );

	forAll(T, celli)
	{
		if (T[celli] == TSat_.value())
		{
			mCon[celli] = 0.0;
			mVap[celli] = 0.0;
		}
		else if (T[celli] < TSat_.value()) mCon[celli] = (T[celli] - TSat_.value())/TSat_.value();
		else							   mVap[celli] = (T[celli] - TSat_.value())/TSat_.value();
	}
	
    return Pair<tmp<volScalarField> >
    (
        (-mcCoeff_)*mCon*(1.0 - limitedAlpha1),

        (-mvCoeff_)*mVap*limitedAlpha1
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Lee::mDotAlphal() const
{
    const volScalarField& T = alpha1_.db().lookupObject<volScalarField>("T");
//    volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));

    volScalarField mCon
    (
        IOobject
        (
            "mCon",
            U_.time().timeName(),
            U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("mCon", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0)
    );
    
    volScalarField mVap
    (
        IOobject
        (
            "mVap",
            U_.time().timeName(),
            U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("mVap", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0)
    );

	forAll(T, celli)
	{
		if (T[celli] == TSat_.value())
		{
			mCon[celli] = 0.0;
			mVap[celli] = 0.0;
		}
		else if (T[celli] < TSat_.value()) mCon[celli] = (T[celli] - TSat_.value())/TSat_.value();
		else							   mVap[celli] = (T[celli] - TSat_.value())/TSat_.value();
	}
	
    return Pair<tmp<volScalarField> >
    (
        (-mcCoeff_)*mCon,

        (-mvCoeff_)*mVap
    );
//    return Pair<tmp<volScalarField> >
//    (
//        (-mcCoeff_)*min((T - TSat_)/TSat_, T0_),
//
//        (-mvCoeff_)*max((T - TSat_)/TSat_, T0_)
//    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Lee::mDotP() const
{
    const volScalarField& T = alpha1_.db().lookupObject<volScalarField>("T");
    volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));

    volScalarField mCon
    (
        IOobject
        (
            "mCon",
            U_.time().timeName(),
            U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("mCon", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0)
    );
    
    volScalarField mVap
    (
        IOobject
        (
            "mVap",
            U_.time().timeName(),
            U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("mVap", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0)
    );

	forAll(T, celli)
	{
		if (T[celli] == TSat_.value())
		{
			mCon[celli] = 0.0;
			mVap[celli] = 0.0;
		}
		else if (T[celli] < TSat_.value()) mCon[celli] = (T[celli] - TSat_.value())/TSat_.value();
		else							   mVap[celli] = (T[celli] - TSat_.value())/TSat_.value();
	}

    return Pair<tmp<volScalarField> >
    (
        (-mcCoeff_)*(1.0 - limitedAlpha1)*mCon,

        (-mvCoeff_)*limitedAlpha1*mVap
    );
//    return Pair<tmp<volScalarField> >
//    (
//        (-mcCoeff_)*(1.0 - limitedAlpha1)*(T - TSat_)/TSat_*neg(T - TSat_),
//
//        (-mvCoeff_)*limitedAlpha1*(T - TSat_)/TSat_*pos(T - TSat_)
//    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Lee::mDotT() const
{
    const volScalarField& T = alpha1_.db().lookupObject<volScalarField>("T");
    volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));

    volScalarField mCon
    (
        IOobject
        (
            "mCon",
            U_.time().timeName(),
            U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("mCon", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0)
    );
    
    volScalarField mVap
    (
        IOobject
        (
            "mVap",
            U_.time().timeName(),
            U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("mVap", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0)
    );

	forAll(T, celli)
	{
		if (T[celli] == TSat_.value())
		{
			mCon[celli] = 0.0;
			mVap[celli] = 0.0;
		}
		else if (T[celli] < TSat_.value()) mCon[celli] = (T[celli] - TSat_.value())/TSat_.value();
		else							   mVap[celli] = limitedAlpha1[celli];
	}

    return Pair<tmp<volScalarField> >
    (
        (-mcCoeff_)*(1.0 - limitedAlpha1)*mCon,

        (-mvCoeff_)*mVap
    );
//    return Pair<tmp<volScalarField> >
//    (
//        (-mcCoeff_)*(1.0 - limitedAlpha1)*(T - TSat_)/TSat_*neg(T - TSat_),
//
//        (-mvCoeff_)*limitedAlpha1*pos(T - TSat_)
//    );
}

void Foam::phaseChangeTwoPhaseMixtures::Lee::correct()
{}


bool Foam::phaseChangeTwoPhaseMixtures::Lee::read()
{
    if (phaseChangeTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");

        phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cc") >> Cc_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cv") >> Cv_;

        mcCoeff_ = Cc_*rho2();
        mvCoeff_ = Cv_*rho1();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
