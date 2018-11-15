    // Solve the momentum equation

    MRF.correctBoundaryVelocity(U);

	tmp<volScalarField> tnut = turbulence->nut();

    // source term M
    tmp<volVectorField> tM
    (
		superFluid->M1()/rho
      + (turbulence->nu() + tnut())*superFluid->M2()
    );
	volVectorField& M = tM.ref();

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
				fvc::interpolate(rhok)*(g & mesh.Sf())
              - fvc::snGrad(p)*mesh.magSf()
              - (fvc::interpolate(M) & mesh.Sf())
            )
        );

        fvOptions.correct(U);
    }
