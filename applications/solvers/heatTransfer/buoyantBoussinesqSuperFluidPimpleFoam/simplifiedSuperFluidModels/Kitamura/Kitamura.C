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

#include "fvc.H"
#include "Kitamura.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace simplifiedSuperFluids
{
    defineTypeNameAndDebug(Kitamura, 0);
    addToRunTimeSelectionTable(simplifiedSuperFluid, Kitamura, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simplifiedSuperFluids::Kitamura::Kitamura
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    simplifiedSuperFluid(typeName, U, phi),
//	KitamuraCoeffs_(subDict(type() + "Coeffs")),
    G_
    (
        IOobject
        (
            "G",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
		U.mesh(),
		dimensionedVector("G", dimensionSet(0,-1,0,1,0,0,0), vector(0,0,0))
    ),

    magG_
    (
        IOobject
        (
            "magG",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
		U.mesh(),
		dimensionedScalar("minDt", dimTemperature/dimLength, SMALL)
    ),

    B_
    (
        IOobject
        (
            "B",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
		U.mesh(),
		dimensionedScalar("B", dimensionSet(0,2,-1,-1,0,0,0), scalar(0))
    )
//    vn_
//    (
//        IOobject
//        (
//            "vn_",
//            U.mesh().time().timeName(),
//            U.mesh(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//		U.mesh(),
//		dimensionedVector("vn_", dimensionSet(0,1,-1,0,0,0,0), vector(0,0,0))
//    ),
//    vs_
//    (
//        IOobject
//        (
//            "vs_",
//            U.mesh().time().timeName(),
//            U.mesh(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//		U.mesh(),
//		dimensionedVector("vs_", dimensionSet(0,1,-1,0,0,0,0), vector(0,0,0))
//    )
{
//    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::simplifiedSuperFluids::Kitamura::calcG()
{
	G_ = fvc::grad(T_);
}

const Foam::volVectorField& Foam::simplifiedSuperFluids::Kitamura::G() const
{
	return G_;
}

const Foam::volScalarField& Foam::simplifiedSuperFluids::Kitamura::B() const
{
	return B_;
}

Foam::tmp<Foam::volVectorField> Foam::simplifiedSuperFluids::Kitamura::calcUn() const
{
	const volVectorField& U = T_.db().lookupObject<volVectorField>("U");
	return tmp<volVectorField>
		   (
		       new volVectorField
			   (
			       IOobject
				   (
				       T_.mesh().time().timeName(),
					   T_.mesh(),
					   IOobject::NO_READ,
					   IOobject::NO_WRITE
				   ),
		           (U - rhos_/rhoHe()*B_*G_)
			   )
		   );
}

Foam::tmp<Foam::volVectorField> Foam::simplifiedSuperFluids::Kitamura::calcUs() const
{
	const volVectorField& U = T_.db().lookupObject<volVectorField>("U");
	return tmp<volVectorField>
		   (
		       new volVectorField
			   (
			       IOobject
				   (
				       T_.mesh().time().timeName(),
					   T_.mesh(),
					   IOobject::NO_READ,
					   IOobject::NO_WRITE
				   ),
		           (U + rhon_/rhoHe()*B_*G_)
			   )
		   );
}

Foam::tmp<Foam::volVectorField> Foam::simplifiedSuperFluids::Kitamura::M1() const
{
	return tmp<volVectorField>
		   (
		       fvc::div(rhon_*rhos_/rhoHe()*B_*B_*G_*G_, "div(V)")
		   );
}

Foam::tmp<Foam::volVectorField> Foam::simplifiedSuperFluids::Kitamura::M2() const
{
	return tmp<volVectorField>
		   (
		       rhos_/rhoHe()*B_*
			   (
			       fvc::laplacian(G_) 
				 + 1./3*fvc::grad(fvc::div(G_))
			   )
    //  (4./3*fvc::laplacian(G) + 1./3*fvc::curl(fvc::curl(G)))
		   );
}


Foam::tmp<Foam::volScalarField> Foam::simplifiedSuperFluids::Kitamura::alpha() const
{
	return tmp<volScalarField> 
		   (
			   pow
		   	   (
		   	       max
		   	       (
		   	           viscosityModelPtr_->k()/magG_/magG_, 
		   	    	   dimensionedScalar("small", dimensionSet(3,3,-9,-3,0,0,0), SMALL)
		   	       ), 
		   	       1./3
		   	   )
		   );
}

Foam::tmp<Foam::volScalarField> Foam::simplifiedSuperFluids::Kitamura::GM() const
{
	return tmp<volScalarField>
		   (
		       new volScalarField
			   (
			       IOobject
				   (
				       T_.mesh().time().timeName(),
					   T_.mesh(),
					   IOobject::NO_READ,
					   IOobject::NO_WRITE
				   ),
				   T_.mesh(),
				   dimensionedScalar("GMzero", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
			   )
		   );
}

Foam::tmp<Foam::volVectorField> Foam::simplifiedSuperFluids::Kitamura::Grest() const
{
	return tmp<volVectorField>
		   (
		       new volVectorField
			   (
			       IOobject
				   (
				       T_.mesh().time().timeName(),
					   T_.mesh(),
					   IOobject::NO_READ,
					   IOobject::NO_WRITE
				   ),
				   T_.mesh(),
		           dimensionedVector("zero", dimensionSet(0,-1,0,1,0,0,0), vector(0,0,0))
			   )
		   );
}

void Foam::simplifiedSuperFluids::Kitamura::correct()
{
	Foam::simplifiedSuperFluid::correct();
	Info<< "Updating G, magG and B..." << endl;
	calcG();
	magG_ = max(mag(G_), dimensionedScalar("minDt", dimTemperature/dimLength, SMALL));
	B_ = pow 
		 (
		     max
		     (
		         sHe()/viscosityModelPtr_->AGM()/rhon_/magG_/magG_, 
		         dimensionedScalar("small", dimensionSet(0,6,-3,-3,0,0,0), SMALL)
		     ), 
		     1./3
		 );
}


bool Foam::simplifiedSuperFluids::Kitamura::read()
{
//    if (simplifiedSuperFluid::read())
//    {
////        KitamuraCoeffs_ = subDict(type() + "Coeffs");
////
////        simplifiedSuperFluidCoeffs_.lookup("Cc") >> Cc_;
////        simplifiedSuperFluidCoeffs_.lookup("Cv") >> Cv_;
////
        return true;
//    }
//    else
//    {
//        return false;
//    }
}


// ************************************************************************* //
