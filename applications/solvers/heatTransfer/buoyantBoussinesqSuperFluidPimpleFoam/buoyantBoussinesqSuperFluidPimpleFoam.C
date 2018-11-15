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

Application
    buoyantBoussinesqSuperFluidPimpleFoam

Group
    grpHeatTransferSolvers

Description
    Transient solver for buoyant, turbulent flow of incompressible 
	superfluid helium.

    Uses the Boussinesq approximation:
    \f[
        rho_{k} = 1 - beta(T - T_{ref})
    \f]

    where:
        \f$ rho_{k} \f$ = the effective (driving) kinematic density
        beta = thermal expansion coefficient [1/K]
        T = temperature [K]
        \f$ T_{ref} \f$ = reference temperature [K]

    Valid when:
    \f[
        \frac{beta(T - T_{ref})}{rho_{ref}} << 1
    \f]

    The solver is based on buoyantBoussinesqPimpleFoam.


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "subCycle.H"
#include "simplifiedSuperFluid.H"
#include "turbulentTransportModel.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "createFields.H"
    #include "createIncompressibleRadiationModel.H"
    #include "createMRF.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "createTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "TControls.H"

            #include "UEqn.H"
            #include "TEqnSubCycle.H"
//            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

		// it is needed for superFluidWall BC
		rhon = superFluid->rhon();
		rhos = superFluid->rhos();
		G = superFluid->G();
		B = superFluid->B();

		// updates q
		s = superFluid->sHe();
		const volScalarField magUn = mag(Un);
		const volScalarField magUs = mag(Us);
		//q = fvc::interpolate(rhos*s*T*(magUn - magUs));
		q = fvc::interpolate(rhos*s*T*(Un.component(1) - Us.component(1)));
		const surfaceScalarField::Boundary& patchHeatFlux = q.boundaryField();
        Info<< "\nWall heat fluxes " << endl;
        forAll(patchHeatFlux, patchi)
        {
           if (typeid(mesh.boundary()[patchi]) == typeid(wallFvPatch))
            {
                Info<< mesh.boundary()[patchi].name()
                    << ": Total "
                    << gSum
                       (
                           mesh.magSf().boundaryField()[patchi]
                          *patchHeatFlux[patchi]
                       )
                    << " [W] over "
                    << gSum
                       (
                           mesh.magSf().boundaryField()[patchi]
                       )
                    << " [m2] ("
                    << gSum
                       (
                           mesh.magSf().boundaryField()[patchi]
                          *patchHeatFlux[patchi]
                       )/
                       gSum 
                       (
                           mesh.magSf().boundaryField()[patchi]
                       )
                    << " [W/m2])"
                    << endl;
            }
	}



        surfaceScalarField gradT=fvc::snGrad(T);
        surfaceScalarField gradp
        (
            IOobject
            (
                "gradp0",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("gradp0", dimTemperature/dimLength, 0)
        );
		if (superFluid->type() == "Pietrowicz")
		{
			gradp = -fvc::snGrad(p*rho)/fvc::interpolate(rho*s);
		}

        surfaceScalarField heatFluxFromAlphaEff =
			fvc::interpolate(alphaEff*cp*rho)*(gradT + gradp);

        const surfaceScalarField::Boundary& patchHeatFlux2 =
                 heatFluxFromAlphaEff.boundaryField();

        Info<< "\nWall heat fluxes from alphaEff" << endl;
        forAll(patchHeatFlux2, patchi)
        {
           if (typeid(mesh.boundary()[patchi]) == typeid(wallFvPatch))
            {
                Info<< mesh.boundary()[patchi].name()
                    << ": Total "
                    << gSum
                       (
                           mesh.magSf().boundaryField()[patchi]
                          *patchHeatFlux2[patchi]
                       )
                    << " [W] over "
                    << gSum
                       (
                           mesh.magSf().boundaryField()[patchi]
                       )
                    << " [m2] ("
                    << gSum
                       (
                           mesh.magSf().boundaryField()[patchi]
                          *patchHeatFlux2[patchi]
                       )/
                       gSum 
                       (
                           mesh.magSf().boundaryField()[patchi]
                       )
                    << " [W/m2])"
                    << endl;
            }
      }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
