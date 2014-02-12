#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <new>
#include "array_structs.h"
#include "numerics.h"
#include "junction_hamiltonian.h"

double H_Kinetic(double pp, double qq, Phys_Parameters &P)
{
	return pp*pp/(2.0*P.m) + qq*qq/(2.0*P.M);
}

double V_M (double ZZ, Phys_Parameters &P)
{
	return P.epsilon * (pow(P.r_equil/ZZ,12) - 2*pow(P.r_equil/ZZ,6));
}

double W_L (double zz, Phys_Parameters &P)
{
	return 1.0/(1+exp(P.alpha*(zz-P.z_L)));
}

double W_R (double zz, Phys_Parameters &P)
{
	return 1.0/(1+exp(-1.0*P.alpha*(zz-P.z_R)));
}

double V_Vac (double zz, Phys_Parameters &P)
{
	double dist = P.z_R - P.z_L;
	return P.e_charge*P.V_b*(1.0 - (zz - P.z_L)/dist)/27.21; //output in atomic units
}

double V_e (double zz, Phys_Parameters &P)
{
	double V_R = P.work_fxn;
	double V_L = P.work_fxn + (P.e_charge * P.V_b)/27.21; //in a.u.
	return W_R(zz,P) * V_R + W_L(zz,P) * V_L + (1.0 - W_R(zz,P) - W_L(zz,P))*V_Vac(zz,P);
}

double V_eMS (double ZZ, Phys_Parameters &P)
{
	if (P.e_affinity == 0)
	{
		return 0.0;
	}
	else
	{
		return P.e_affinity - (P.cconst/(27.21*(ZZ-P.z_L)));
	}
}

double W_coup (double zz, double ZZ, Phys_Parameters &P)
{
	if (P.coupling == 0)
	{
		return exp(-1.0*(zz-ZZ)*(zz-ZZ)/(2.0*P.sigma*P.sigma));
	}
	else
	{
		double diff 	= abs(zz-ZZ);
		double parab 	= (diff/P.rad) * (diff/P.rad) - 1.0;
		double on 		= 1.0/(1.0 + exp(P.del*(diff-P.rad)));
		double off 		= 1.0/(1.0 + exp(-1.0*P.del*(diff+P.rad)));
		return -1.0*parab*on*off;
	}
}

double Elec_parab (double zz, double ZZ, Phys_Parameters &P)
{
	double diff = abs(zz-ZZ);
	return (diff/P.rad)*(diff/P.rad)*abs(P.e_affinity);
}

double V_eM (double zz, double ZZ, Phys_Parameters &P)
{
	return W_coup(zz,ZZ,P)*(V_eMS(ZZ,P) - V_e(zz,P));
}

double V_0 (double zz, double ZZ, Phys_Parameters &P)
{
	return V_M(ZZ,P) + V_e(zz,P);
}

double absorbing_potential(double XX, Phys_Parameters &P, Comp_Parameters &C)
{
	double c = 2.622;
	double xx;
	if (XX <= (C.xmin + P.Lbarrier))
	{
		xx = -1.0 * (c +0.001)*(XX-(C.xmin + P.Lbarrier))/P.Lbarrier;
		return P.Lbarrier_itn*(4.0/((c-xx)*(c-xx))+4.0/((c+xx)*(c+xx))-8.0/(c*c))*(cosh(log(2.0+sqrt(3.0))*xx/c)-1.0);
	}
	else if (XX >= (C.xmax - P.Rbarrier))
	{
		xx = (c+0.001)*(XX-(C.xmax - P.Rbarrier))/P.Rbarrier;
		return P.Rbarrier_itn*(4.0/((c-xx)*(c-xx))+4.0/((c+xx)*(c+xx))-8.0/(c*c))*(cosh(log(2.0+sqrt(3.0))*xx/c)-1.0);
	}
	else
	{
		return 0;
	}
}

//	Absorbing potential in the molecular dimension
double absorbing_potential_mol(double XX, Comp_Parameters &C)
{
	double c = 2.622;
	double xx;
	double Lbarrier = 3.0;
	double Rbarrier = 5.0;
	if (XX <= (C.ymin + Lbarrier))
	{
		xx = -1.0 * (c +0.001)*(XX-(C.ymin + Lbarrier))/Lbarrier;
		return (2*M_PI*M_PI/(Lbarrier*Lbarrier))*(4.0/((c-xx)*(c-xx))+4.0/((c+xx)*(c+xx))-8.0/(c*c))*(cosh(log(2.0+sqrt(3.0))*xx/c)-1.0);
	}
	else if (XX >= (C.ymax - Rbarrier))
	{
		xx = (c+0.001)*(XX-(C.ymax - Rbarrier))/Rbarrier;
		return (2*M_PI*M_PI/(Rbarrier*Rbarrier))*(4.0/((c-xx)*(c-xx))+4.0/((c+xx)*(c+xx))-8.0/(c*c))*(cosh(log(2.0+sqrt(3.0))*xx/c)-1.0);
	}
	else
	{
		return 0;
	}
}


double H_Potential_NoAbs (double zz, double ZZ, Phys_Parameters &P)
{
	return (V_M(ZZ,P) + V_e(zz,P) + V_eM(zz,ZZ,P));
}

cplx H_Potential (double zz, double ZZ, Phys_Parameters &P, Comp_Parameters &C)
{
	return -1.0j*(absorbing_potential(zz,P,C) + absorbing_potential_mol(ZZ,C))+ H_Potential_NoAbs(zz,ZZ,P);
}

double wvfxn_nuclear (double ZZ, Phys_Parameters &P)
{
	return ho_wvfxn(0, ZZ, P.r_equil,P.omega,P.M,P.hbar);
}

cplx wvfxn_electron (double zz, Phys_Parameters &P)
{
	double coeff 	= sqrt(sqrt(1.0/(M_PI*P.xi*P.xi)));
	cplx arg1 		= 1.0j * P.kvect * zz / P.hbar;
	arg1 			= exp (arg1);
	double arg2 	= exp (-1.0 * (zz-P.z0) * (zz-P.z0)/(2.0*P.xi*P.xi));
	cplx val 		= coeff * arg1 * arg2;
	return val;
}
