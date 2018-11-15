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

#include "simplifiedSuperFluid.H"
#include "Helium.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simplifiedSuperFluid, 0);
    defineRunTimeSelectionTable(simplifiedSuperFluid, components);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simplifiedSuperFluid::simplifiedSuperFluid
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    singlePhaseTransportModel(U, phi),
    simplifiedSuperFluidCoeffs_(subDict(type + "Coeffs")),
	T_(U.db().lookupObject<volScalarField>("T")),
	//rhoHe_(U.db().lookupObject<volScalarField>("rho")),
    //rhoHe_
    //(
    //    IOobject
    //    (
    //        "rhoHe",
    //        U.mesh().time().timeName(),
    //        U.mesh(),
    //        IOobject::NO_READ,
    //        IOobject::NO_WRITE
    //    ),
    //    U.mesh(),
	//	dimensionedScalar("rhoHe", dimensionSet(1,-3,0,0,0,0,0), 0.0)
    //),
    rhon_
    (
        IOobject
        (
            "rhon",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
		dimensionedScalar("rhon0", dimensionSet(1,-3,0,0,0,0,0), 0.0)
    ),

    rhos_
    (
        IOobject
        (
            "rhos",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
		dimensionedScalar("rhos", dimensionSet(1,-3,0,0,0,0,0), 0.0)
    ),

    Pr_
    (
        IOobject
        (
            "Pr",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
    	dimensionedScalar("Pr0", dimless, 0.0)
    )
{

}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::simplifiedSuperFluid::correct()
{
	viscosityModelPtr_->correct();
	rhons();
	calcPr();
}

void Foam::simplifiedSuperFluid::rhons()
{
	const dimensionedScalar Tlambda = Foam::viscosityModels::Helium::Tlambda();
	rhon_ = rhoHe()*pow(max(T_/Tlambda, dimensionedScalar("small", dimless, SMALL)), scalar(5.6));
	rhos_ = rhoHe() - rhon_;
}

bool Foam::simplifiedSuperFluid::read()
{
    if (singlePhaseTransportModel::read())
    {
        simplifiedSuperFluidCoeffs_ = subDict(type() + "Coeffs");
//        lookup("pSat") >> pSat_;
//        lookup("cutoff") >> cutoff_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
