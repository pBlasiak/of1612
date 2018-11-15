{
    surfaceScalarField alpha1f =
        min(max(fvc::interpolate(alpha1), scalar(0)), scalar(1));

    surfaceScalarField kf =
            alpha1f*k1 + (scalar(1) - alpha1f)*k2;
    surfaceScalarField cpf =
            alpha1f*cp1 + (scalar(1) - alpha1f)*cp2;
	volScalarField cp = alpha1*cp1 + (1.0 - alpha1)*cp2;


    volScalarField IF
    (
        IOobject
        (
            "IF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
		mesh,
        dimensionedScalar("IF", dimensionSet(1,-1,-3,-1,0,0,0), scalar(0))
    );

    Pair<tmp<volScalarField> > mDotT =
        mixture->mDotT();
    const volScalarField& mDotTc = mDotT[0]();
    const volScalarField& mDotTv = mDotT[1]();

//	IF = twoPhaseProperties->IF();
	forAll(alpha1, celli)
	{
		if (alpha1[celli] < 1.0 - cutoff.value() && alpha1[celli] > cutoff.value()) IF[celli] = 1.0;
	}

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(rhoPhi, T)
      - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), T)
      - fvm::laplacian(kf/cpf, T)
     ==
	    hEvap/cp*(fvm::Sp(mDotTv/TSat, T) + mDotTc - mDotTv)  // tak jest jak m = mc + mv
//	    hEvap/cp*(mDotTc - mDotTv*(T - TSat)/TSat)  // tak jest jak m = mc - mv
//	  - 1./cp*fvm::Sp(pow(10,10)*IF, T) + 1./cp*pow(10,10)*TSat*IF
    );


//    TEqn.relax();
    TEqn.solve();

    Info<< "min/max(T) = " << min(T).value() << ", "
        << max(T).value() <<endl;
}
