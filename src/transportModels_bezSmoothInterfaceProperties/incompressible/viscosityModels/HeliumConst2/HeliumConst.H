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

Class
    Foam::viscosityModels::HeliumConst

Description
     Standard power-law non-Newtonian viscosity model.

SourceFiles
    HeliumConst.C

\*---------------------------------------------------------------------------*/

#ifndef HeliumConst_H
#define HeliumConst_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class HeliumConst Declaration
\*---------------------------------------------------------------------------*/

class HeliumConst
:
    public viscosityModel
{
    // Private data

		//- lambda Temperature 2.1711132461
		static const dimensionedScalar Tlambda_;

        dictionary HeliumConstCoeffs_;

        //- Thermal expansion coefficient
        dimensionedScalar betaHe0_;

        //- Thermal expansion coefficient
        volScalarField betaHe_;
		
        //- Goerter-Mellink parameter
        dimensionedScalar AGMHe0_;

        //- Goerter-Mellink parameter
        volScalarField AGMHe_;

        //- Entropy
        dimensionedScalar sHe0_;
		
        //- Entropy
        volScalarField sHe_;
		
        //- Specific heat at constant pressure
        dimensionedScalar cpHe0_;

        //- Specific heat at constant pressure
        volScalarField cpHe_;
		
        //- Thermal conductivity function
        dimensionedScalar onebyf0_;

        //- Thermal conductivity function
        volScalarField onebyf_;

        //- Density of the mixture 
        dimensionedScalar rhoHe_;

        dimensionedScalar nu0_;

        volScalarField nu_;



    // Private Member Functions

public:

    //- Runtime type information
    TypeName("HeliumConst");


    // Constructors

        //- Construct from components
        HeliumConst
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    ~HeliumConst()
    {}


    // Member Functions

        //- Return the laminar viscosity
        tmp<volScalarField> nu() const override
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const override
        {
            return nu_.boundaryField()[patchi];
        }

		//- Returns thermal expansion coefficient of He
        const volScalarField& beta() const override
		{
			return betaHe_;
		}

		//- Returns Goerter-Mellink coefficient for helium based on temperature
        const volScalarField& AGM() const override
		{
			return AGMHe_;
		}

		//- Returns entropy for helium based on temperature
        const volScalarField& s() const override
		{
			return sHe_;
		}

		//- Returns dynamic viscosity for helium based on temperature
        const volScalarField& eta() const override
		{
			NotImplemented;
			return nu_*rhoHe_;
		}

		//- Returns specific heat capacity for helium based on temperature
        const volScalarField& cp() const override
		{
			return cpHe_;
		}

		//- Returns thermal conductivity function for helium based on temperature
        const volScalarField& onebyf() const override
		{
			return onebyf_;
		}

		const dimensionedScalar& rho() const override
		{
			return rhoHe_;
		}

		static const dimensionedScalar Tlambda() 
		{
			return Tlambda_;
		}

        //- Correct the laminar viscosity and thermal properties of He
        void correct()
		{}

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
