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

#include "superFluidWallUsFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::superFluidWallUsFvPatchVectorField::
superFluidWallUsFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::superFluidWallUsFvPatchVectorField::
superFluidWallUsFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::superFluidWallUsFvPatchVectorField::
superFluidWallUsFvPatchVectorField
(
    const superFluidWallUsFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::superFluidWallUsFvPatchVectorField::
superFluidWallUsFvPatchVectorField
(
    const superFluidWallUsFvPatchVectorField& sfwpvf
)
:
    fixedValueFvPatchVectorField(sfwpvf)
{}


Foam::superFluidWallUsFvPatchVectorField::
superFluidWallUsFvPatchVectorField
(
    const superFluidWallUsFvPatchVectorField& sfwpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(sfwpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::superFluidWallUsFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	// Accessing face normal vectors
    const vectorField nf(patch().nf());

	// Accessing G vector at the patch
	const fvPatchField<vector>& Gwall = patch().lookupPatchField<volVectorField, vector>("G");

    // Calculated the component of the G vector normal to the wall
    vectorField GwallNorm = (nf & Gwall)*nf;
	
    // Calculated the component of the G vector parallel to the wall
    vectorField GwallTang = Gwall - (nf & Gwall)*nf;
	//Info<< patch().name() << endl;
	//Info<< GwallNorm << endl;

	// Accessing other variables at the patch
	const fvPatchField<scalar>& rhonWall = patch().lookupPatchField<volScalarField, scalar>("rhon");
	const fvPatchField<scalar>& Bwall = patch().lookupPatchField<volScalarField, scalar>("B");
	const fvPatchField<scalar>& rhoWall = patch().lookupPatchField<volScalarField, scalar>("rho");

    vectorField UsNorm = rhonWall/rhoWall*Bwall*GwallNorm;
    vectorField UsTang = Bwall*GwallTang;
    vectorField UsTotal = UsNorm + UsTang;

    vectorField::operator=(UsTotal);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::superFluidWallUsFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        superFluidWallUsFvPatchVectorField
    );
}

// ************************************************************************* //
