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
#include "HeliumConst.H"

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
    public HeliumConst
{
    // Protected data
protected:

        dictionary HeliumCoeffs_;

        const volScalarField& T_;

        dimensionedScalar nuMin_;
        dimensionedScalar nuMax_;




    // Private Member Functions

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

        //- Correct the laminar viscosity and thermal properties of He
        void correct() override;

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
