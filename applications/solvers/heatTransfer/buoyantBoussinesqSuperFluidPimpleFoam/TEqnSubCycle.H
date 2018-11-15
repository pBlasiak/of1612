{
    if (nTSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
        volScalarField rhokSum("rhokSum", rhok-rhok);

        for
        (
            subCycle<volScalarField> TSubCycle(T, nTSubCycles);
            !(++TSubCycle).end();
        )
        {
            #include "TEqn.H"
            rhokSum += (runTime.deltaT()/totalDeltaT)*rhok;
        }

        rhok = rhokSum;
    }
    else
    {
        #include "TEqn.H"
    }

//    rho == alpha1*rho1 + alpha2*rho2;
}
