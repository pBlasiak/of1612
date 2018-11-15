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
    Foam::viscosityModels::Helium

Description
     Standard power-law non-Newtonian viscosity model.

SourceFiles
    Helium.C

\*---------------------------------------------------------------------------*/

#ifndef Helium_H
#define Helium_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class Helium Declaration
\*---------------------------------------------------------------------------*/

class Helium
:
    public viscosityModel
{
    // Private data

		//- lambda Temperature 2.1711132461
		static const dimensionedScalar Tlambda_;

		//- Minimum temperature TMin = 1.5 K
		static const dimensionedScalar TMin_;

		//- Maximum temperature TMax = 2.167 K
		static const dimensionedScalar TMax_;

		//- Minimum index number for thermophysical tables i = 0
		static const label indexMin_;

		//- Maximum index number for thermophysical tables i = 667
		static const label indexMax_;

		//- Step of temperature in the thermophysical tables dT = 0.001 K
		static const scalar dT_;

		//- Table for betaHe
		static const List<scalar> betaHeTable_;

		//- Table for AGMHe
		static const List<scalar> AGMHeTable_;

		//- Table for sHe
		static const List<scalar> sHeTable_;

		//- Table for etaHe
		static const List<scalar> etaHeTable_;

		//- Table for cpHe
		static const List<scalar> cpHeTable_;

		//- Table for onebyf
		static const List<scalar> onebyfTable_;

		//- Table for rhoHe
		static const List<scalar> rhoHeTable_;


        dictionary HeliumCoeffs_;

        //- Thermal expansion coefficient
        volScalarField betaHe_;
		
        //- Goerter-Mellink parameter
        volScalarField AGMHe_;

        //- Entropy
        volScalarField sHe_;
		
        //- Dynamic viscosity
        volScalarField etaHe_;
		
        //- Specific heat at constant pressure
        volScalarField cpHe_;
		
        //- Thermal conductivity function
        volScalarField onebyf_;

        //- Density of the mixture 
        dimensionedScalar rhoHe_;

        const volScalarField& T_;

        dimensionedScalar nuMin_;
        dimensionedScalar nuMax_;

        volScalarField nu_;



    // Private Member Functions

//		//- Calculates helium dynamic viscosity based on temperature
//        void calcEtaHe();

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() ;

		//- Calculates helium thermal property based on temperature
        void calcHeProp(volScalarField& vsf, const List<scalar>& vsfTable);


public:

    //- Runtime type information
    TypeName("Helium");


    // Constructors

        //- Construct from components
        Helium
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    ~Helium()
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
			return etaHe_;
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
        void correct();

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
