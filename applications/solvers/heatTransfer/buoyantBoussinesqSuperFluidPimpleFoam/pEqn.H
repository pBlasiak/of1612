{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

    surfaceScalarField phig(rAUf*fvc::interpolate(rhok)*(g & mesh.Sf()));
    surfaceScalarField phiM(-rAUf*(fvc::interpolate(M) & mesh.Sf()));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + rAUf*fvc::ddtCorr(U, phi)
//      + phig
//      + phiM
    );

    MRF.makeRelative(phiHbyA);
	adjustPhi(phiHbyA, U, p);
	phiHbyA += phig + phiM;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAUf, MRF);

   // nie wiem czy to jest potrzebne
   // phiHbyA -= phiM;                                                     //3  (+)

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, getRefCellValue(p, pRefCell));

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi = phiHbyA - pEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U = HbyA + rAU*fvc::reconstruct((phiM + phig - pEqn.flux())/rAUf); //5   (+)
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

	Un = superFluid->calcUn()();
	Un.correctBoundaryConditions();
	Us = superFluid->calcUs()();
	Us.correctBoundaryConditions();
	UnmUs = Un - Us;
	
    Info<< "min/max(U) = " << min(mag(U)).value() << ", "
        << max(mag(U)).value() <<endl;

//    p = p_rgh + rhok*gh;
//
//    if (p_rgh.needReference())
//    {
//        p += dimensionedScalar
//        (
//            "p",
//            p.dimensions(),
//            pRefValue - getRefCellValue(p, pRefCell)
//        );
//        p_rgh = p - rhok*gh;
//    }
}
