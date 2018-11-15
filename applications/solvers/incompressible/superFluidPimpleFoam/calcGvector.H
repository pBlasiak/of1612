// Calculate source terms

Info<< "Updating G vector, rhon, rhos and B" << endl;

volVectorField gradp(-fvc::grad(p*rho)/rho/s);
volVectorField gradT(fvc::grad(T));
G = gradp + gradT;
//G = gradT;
magG = max(mag(G), dimensionedScalar("minDt", dimTemperature/dimLength, SMALL));

rhon = rho*pow(max(T/Tlambda, dimensionedScalar("small", dimless, SMALL)), scalar(5.6));
rhos = rho - rhon;

B = pow(max(s/AGM/rhon/magG/magG, dimensionedScalar("small", dimensionSet(0,6,-3,-3,0,0,0), SMALL)), 1./3);
