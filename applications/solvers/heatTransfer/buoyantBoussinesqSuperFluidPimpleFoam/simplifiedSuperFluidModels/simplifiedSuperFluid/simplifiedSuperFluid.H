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

Class
    Foam::simplifiedSuperFluid

Description
	Base class for simplified superfluid helium models. 

SourceFiles
    simplifiedSuperFluid.C
    newSimplifiedSuperFluid.C

\*---------------------------------------------------------------------------*/

#ifndef simplifiedSuperFluid_H
#define simplifiedSuperFluid_H

#include "singlePhaseTransportModel.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class simplifiedSuperFluid Declaration
\*---------------------------------------------------------------------------*/

class simplifiedSuperFluid
:
    public singlePhaseTransportModel
{

protected:

    // Protected data

        dictionary simplifiedSuperFluidCoeffs_;

		const volScalarField& T_;

		//const volScalarField& rhoHe_;

        //- Normal density component
        volScalarField rhon_;

        //- Superfluid density component
        volScalarField rhos_;

        //- Superfluid Prandtl number
        volScalarField Pr_;

    // Private Member Functions

        //- Disallow copy construct
        simplifiedSuperFluid(const simplifiedSuperFluid&);

        //- Disallow default bitwise assignment
        void operator=(const simplifiedSuperFluid&);

		//- Updates rhon_ and rhos_
		void rhons();

public:


    //- Runtime type information
    TypeName("simplifiedSuperFluid");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            simplifiedSuperFluid,
            components,
            (
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (U, phi)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<simplifiedSuperFluid> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Constructors

        //- Construct from components
        simplifiedSuperFluid
        (
            const word& type,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~simplifiedSuperFluid()
    {}


    // Member Functions

		//- Updates G vector
		virtual void calcG() = 0;

		//- Calculates normal velocity component un
		virtual tmp<volVectorField> calcUn() const = 0;
	
		//- Calculates super velocity component us
		virtual tmp<volVectorField> calcUs() const = 0;

        //- Returns reference to G vector
        virtual const volVectorField& G() const = 0;

        //- Returns B 
		//  B = (s/A/rhon/magG^2)^(1./3)
        virtual const volScalarField& B() const = 0;

        //- Returns convective momentum equation source term
		// div(rhon*rhos/rho/rho*B*B*G*G);
        virtual tmp<volVectorField> M1() const = 0;

        //- Returns viscous momentum equation source term
		//  It has to be multiplied by nut()
        virtual tmp<volVectorField> M2() const = 0;

        //- Returns laminar thermal conductivity 
		//  alpha = (1/f/|G|^2)^1/3
        virtual tmp<volScalarField> alpha() const = 0;

        //- Returns Goerter-Mellink energy equation source terms
        virtual tmp<volScalarField> GM() const = 0;

        //- Returns the rest part of G vector energy equation source terms
        //virtual tmp<volScalarField> Grest() const = 0;
        virtual tmp<volVectorField> Grest() const = 0;

        //- Correct thermophysical properties of He
        virtual void correct();

        //- Read the transportProperties dictionary and update
        virtual bool read();

		const volScalarField& betaHe() const
		{
			return viscosityModelPtr_->beta();
		}

		const volScalarField& cpHe() const
		{
			return viscosityModelPtr_->cp();
		}

		const volScalarField& sHe() const
		{
			return viscosityModelPtr_->s();
		}

		const volScalarField& etaHe() const
		{
			return viscosityModelPtr_->eta();
		}

		const volScalarField& rhoHe() const
		{
			return viscosityModelPtr_->rho();
		}

		void calcPr()
		{
			Pr_ = viscosityModelPtr_->cp()*viscosityModelPtr_->eta()/alpha();
		}

		const volScalarField& Pr() const
		{
			return Pr_;
		}
		
		const volScalarField& rhos() const
		{
			return rhos_;
		}

		const volScalarField& rhon() const
		{
			return rhon_;
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
