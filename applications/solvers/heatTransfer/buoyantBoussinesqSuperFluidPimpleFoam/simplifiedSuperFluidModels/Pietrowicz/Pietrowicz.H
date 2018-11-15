/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::simplifiedSuperFluids::Pietrowicz

Description
    Pietrowicz cavitation model slightly modified so that the condensation term
    is switched off when the pressure is less than the saturation vapour
    pressure.  This change allows the condensation term to be formulated as
    a coefficient multiplying (p - p_sat) so that it can be included as an
    implicit term in the pressure equation.

    Reference:
    \verbatim
        Pietrowicz, R.F., Boger, D.A., Stinebring, D.R., Chyczewski, Lindau. J.W.,
        Gibeling, H.J., Venkateswaran, S., Govindan, T.R.,
        "A Preconditioned Implicit Method for Two-Phase Flows with Application
         to Cavitation Prediction,"
        Computers and Fluids,
        29(8):849-875, 2000.
    \endverbatim

SourceFiles
    Pietrowicz.C

\*--------------------------------------------------------------------*/

#ifndef Pietrowicz_H
#define Pietrowicz_H

#include "Kitamura.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace simplifiedSuperFluids
{

/*--------------------------------------------------------------------*\
                              Class Pietrowicz
\*--------------------------------------------------------------------*/

class Pietrowicz
:
    public Kitamura
{
    // Private data

		//- Reference to pressure field
		
        dictionary PietrowiczCoeffs_;

        const volScalarField& p_;

		volVectorField Grest_;

		Switch GM_;
//		Switch GRest_;


public:

    //- Runtime type information
    TypeName("Pietrowicz");


    // Constructors

        //- construct from components
        Pietrowicz
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~Pietrowicz()
    {}


    // Member Functions
	
        //- Read the transportProperties dictionary and update
        virtual bool read() override;

		//- Updates G vector
		virtual void calcG() override;

        //- Returns Goerter-Mellink energy equation source terms
        virtual tmp<volScalarField> GM() const override;

        //- Returns the rest part of G vector energy equation source terms
        //virtual tmp<volScalarField> Grest() const override;
        virtual tmp<volVectorField> Grest() const override;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace simplifiedSuperFluids
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
