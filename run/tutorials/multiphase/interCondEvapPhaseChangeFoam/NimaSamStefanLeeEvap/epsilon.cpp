#include <iostream>
#include <cmath>
//#include <cstdlib>
#include <fstream>
#include <cstdio>
/*
Program solves transcendental equation:
	epsilon*exp(epsilon^2)*erf(epsilon) - constVal = 0
with use of Newton's method. Than it calculates deltaInter
which is the position of the interface as a exact solution
of one-dimensional Stefan problem for evaporation.
*/

using namespace std;

double func(double x, const double constValue)
{
	return x*exp(pow(x,2))*erf(x) - constValue;
}

double derivErf(const double x)
{
	return 2*exp(-pow(x,2))/sqrt(M_PI);
}

double derivFunc(double x)
{
	return exp(pow(x,2))*erf(x) + x*(2*x*exp(pow(x,2))*erf(x) + exp(pow(x,2))*derivErf(x) ); 
}

double deltaInter(const double eps, const double thermCond, const double rho, const double cp, const double time)
{
	return 2*eps*sqrt(thermCond*time/rho/cp);
}

int main()
{ 
    std::cout << "\nUsuwam plik interfacePosStefProbl_analSol.dat ...\n" << std::endl;
	remove("interfacePosStefProbl_analSol.dat");

	const long maxIterNr = 20000;
	long iter = 0;
	const double Tw = 398.15;
	const double TSat = 373.15;
	const double cpv = 1007;
	const double rhov = 1;
	const double kappav = 0.0257;
	const double hfg = 2200e+3;
	const double tol = 1e-12;
	double epsilon = 0;
	double epsilonPrev = 0;
	const double LHS = cpv*(Tw - TSat)/hfg/sqrt(M_PI);

	std::cout << "LHS = " << LHS << endl;
	std::cout << "Type the starting value for epsilon:" << std::endl;	
	std::cin >> epsilon;

	do 
	{
		epsilonPrev = epsilon;		
		epsilon -= func(epsilon, LHS)/derivFunc(epsilon);	
		iter++;
	} while ( abs(epsilon - epsilonPrev) > tol && iter < maxIterNr );

    std::cout << "Liczba iteracji: " << iter << endl;
	std::cout << "epsilon = " << epsilon << endl;
	std::cout << "deltaInter(t = 1s) = " << deltaInter(epsilon,kappav, rhov, cpv, 1.) << std::endl;

	double totTime;
	std::cout << "Podaj czas koncowy (liczba dodatnia): " << std::endl;
	std::cin >> totTime; 

	std::fstream plik("interfacePosStefProbl_analSol.dat", std::ios_base::app | std::ios_base::out);
	if( plik.good() == true )
	{
		std::cout << "Podaj krok czasowy: " << std::endl;
		double dt = 1.0;
		std::cin >> dt;
	    std::cout << "\nZapisuje wyniki do pliku interfacePosStefProbl_analSol.dat ...\n" << std::endl;
    	for (double t = 0; t <= totTime; t+=dt)
    	{
    		std::cout << "deltaInter(t = " << t << " s) = " 
				      << deltaInter(epsilon,kappav, rhov, cpv, t) << " m" << std::endl;
    	    plik << t << "\t" << deltaInter(epsilon,kappav, rhov, cpv, t) << std::endl;
    	}
		plik.close();
	}
	else 
	{
	    std::cout << "ERORR! Nie mozna otworzyc pliku!" << std::endl;
	}

	return 0;
}
