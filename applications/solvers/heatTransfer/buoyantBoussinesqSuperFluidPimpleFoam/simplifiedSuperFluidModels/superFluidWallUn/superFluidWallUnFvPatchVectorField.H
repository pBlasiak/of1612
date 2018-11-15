/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::superFluidWallUnFvPatchVectorField

Group
    grpWallBoundaryConditions

Description
    This boundary condition provides a boundary condition for 
	normal component of velocity for cases with
    simplified super fluid helium models. 

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            superFluidWallUn;
        value           uniform (0 0 0);    // Initial value
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchVectorField

SourceFiles
    superFluidWallUnFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef superFluidWallUnFvPatchVectorField_H
#define superFluidWallUnFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class superFluidWallUnFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class superFluidWallUnFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("superFluidWallUn");


    // Constructors

        //- Construct from patch and internal field
        superFluidWallUnFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        superFluidWallUnFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given superFluidWallUnFvPatchVectorField
        //  onto a new patch
        superFluidWallUnFvPatchVectorField
        (
            const superFluidWallUnFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        superFluidWallUnFvPatchVectorField
        (
            const superFluidWallUnFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new superFluidWallUnFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        superFluidWallUnFvPatchVectorField
        (
            const superFluidWallUnFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new superFluidWallUnFvPatchVectorField(*this, iF)
            );
        }


    // Member functions
	
        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
