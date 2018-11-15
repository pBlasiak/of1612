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
    Foam::viscosityModels::HeliumFineOnebyfConstRho

Description
     Standard power-law non-Newtonian viscosity model.

SourceFiles
    HeliumFineOnebyfConstRho.C

\*---------------------------------------------------------------------------*/

#ifndef HeliumFineOnebyfConstRho_H
#define HeliumFineOnebyfConstRho_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "HeliumConstRho.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class HeliumFineOnebyfConstRho Declaration
\*---------------------------------------------------------------------------*/

class HeliumFineOnebyfConstRho
:
    public HeliumConstRho
{
    // Private data

		static const label largeIndex_;

		//- Step of temperature in the thermophysical tables dT = 0.0001 K
		static const scalar dTFineOnebyf_;

		//- Table for betaHe
		static const List<scalar> fineOnebyfTable_;


    // Private Member Functions

		//- Calculates helium 1/f function based on temperature
        void fineOnebyf();


public:

    //- Runtime type information
    TypeName("HeliumFineOnebyfConstRho");


    // Constructors

        //- Construct from components
        HeliumFineOnebyfConstRho
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~HeliumFineOnebyfConstRho()
    {}


    // Member Functions

        //- Correct the laminar viscosity and thermal properties of He
        virtual void correct() override;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
