#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <new>
#include <fftw3.h>
#include "array_structs.h"
#include "numerics.h"

double hermite(int n, double x)
{
	if (n == 0) 
	{
		return 1;
	}
	else if (n == 1) 
	{
		return 2*x;
	}
	else 
	{
		double poly = 2*x*hermite(n-1,x) - 2*(n-1)*hermite(n-2,x);
		return poly;
	}
}

double factorial(int n)
{
	double fact;
	if (n==0) 
	{
		fact = 1;
	}
	else 
	{
		fact = n*factorial(n-1);
	}
	return fact;
}

double ho_wvfxn(int n, double x, double xcen, double omega, double m, double hbar)
{ //nth harmonic oscillator solution evaluated at x, centered at xcen
	double coeff1 	= 1/sqrt(pow(2,n)*factorial(n));
	double coeff2 	= sqrt(sqrt(m*omega/M_PI/hbar));
	double coeff 	= coeff1 * coeff2;
	x 				-= xcen;
	double fxn 		= coeff * exp(-1*m*omega*x*x/2/hbar) * hermite(n,sqrt(m*omega/hbar)*x);
	return fxn;
}

int kron_delta(int x, int y)
{
	if (x == y) 
	{
		return 1;
	} 
	else 
	{
		return 0;
	}
}
