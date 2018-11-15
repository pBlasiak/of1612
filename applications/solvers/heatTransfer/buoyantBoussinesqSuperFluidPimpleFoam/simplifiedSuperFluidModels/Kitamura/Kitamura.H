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
    Foam::simplifiedSuperFluids::Kitamura

Description
    Kitamura super fluid helium model.

    Reference:
    \verbatim
        Kitamura T., Shiramizu K., Fujimoto N., Rao Y.F., Fukuda K.,
        "A numerical model on transient, two-dimensional flow and heat
		transfer in He II", Cryogenics 37:1-9, 1997.
    \endverbatim

SourceFiles
    Kitamura.C

\*--------------------------------------------------------------------*/

#ifndef Kitamura_H
#define Kitamura_H

#include "simplifiedSuperFluid.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace simplifiedSuperFluids
{

/*--------------------------------------------------------------------*\
                              Class Kitamura
\*--------------------------------------------------------------------*/

class Kitamura
:
    public simplifiedSuperFluid
{
    // Protected data
protected:
		
//		dictionary KitamuraCoeffs_;

        volVectorField G_;

		volScalarField magG_;

		//- B = (s/A/rhon/magG^2)^(1./3)
		volScalarField B_;

//        //- Normal velocity component
//        volVectorField vn_;
//
//        //- Super velocity component
//        volVectorField vs_;

public:

    //- Runtime type information
    TypeName("Kitamura");


    // Constructors

        //- construct from components
        Kitamura
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~Kitamura()
    {}


    // Member Functions
	
		//- Updates G vector
		virtual void calcG() override;

		//- Calculates normal velocity component un
		virtual tmp<volVectorField> calcUn() const override;
	
		//- Calculates super velocity component us
		virtual tmp<volVectorField> calcUs() const override;

        //- Returns reference to G vector
        virtual const volVectorField& G() const override;

        //- Returns B 
		//  B = (s/A/rhon/magG^2)^(1./3)
        virtual const volScalarField& B() const override;

        //- Returns convective momentum equation source term
		// div(rhon*rhos/rho/rho*B*B*G*G);
        virtual tmp<volVectorField> M1() const override;

        //- Returns viscous momentum equation source term
		//  It has to be multiplied by nut()
        virtual tmp<volVectorField> M2() const override;

        //- Returns laminar thermal conductivity 
		//  alpha = (1/f/|G|^2)^1/3
        virtual tmp<volScalarField> alpha() const override;

        //- Returns Goerter-Mellink energy equation source terms
        virtual tmp<volScalarField> GM() const override;

        //- Returns the rest part of G vector energy equation source terms
       // virtual tmp<volScalarField> Grest() const override;
        virtual tmp<volVectorField> Grest() const override;

        //- Correct variables G, magG, B
        virtual void correct() override;

        //- Read the transportProperties dictionary and update
        virtual bool read() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace simplifiedSuperFluids
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
