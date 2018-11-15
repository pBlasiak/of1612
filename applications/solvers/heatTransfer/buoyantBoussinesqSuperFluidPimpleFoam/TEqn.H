{
    alphat = tnut()/Prt;
    alphat.correctBoundaryConditions();

    tmp<volScalarField> talphaEff
    (
	    (superFluid->alpha()/rho/cp + alphat)
    );
	alphaEff = talphaEff.ref();


	tmp<volScalarField> tGMpGrest
	(
		superFluid->GM()/rho/cp + fvc::div(alphaEff*superFluid->Grest(), "div(Grest)")
	);
	volScalarField& GMpGrest = tGMpGrest.ref();

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
//      - fvm::Sp(fvc::div(phi),T) // Taken from Nima Sam (post nr 107)
     ==
        radiation->ST(rhoCpRef, T) // nie policzone rhoCpRef
	  + GMpGrest
//      + fvm::SuSp(GMpGrest/T, T)
      + fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    radiation->correct();

    fvOptions.correct(T);
    T.correctBoundaryConditions();

	superFluid->correct();

	beta = superFluid->betaHe();
    rhok = 1.0 - beta*(T - TRef);
	cp = superFluid->cpHe();
	rho = superFluid->rhoHe();
	Pr = superFluid->Pr(); // it is needed for alphat BC

    Info<< "min/max(T) = " << min(T).value() << ", "
        << max(T).value() <<endl;
}
