#include "junction.hpp"

using namespace std;

double HKinetic1D(double pp, double mass)
{
  return pp*pp/(2.0*mass);
} ///Kinetic energy operator for one dimension

double HKinetic2D(double pp, double qq, double mass1, double mass2)
{
  return pp*pp/(2.0*mass1) + qq*qq/(2.0*mass2);
}

double VM (double ZZ, programInputs &IP)
{
  return IP.epsilon * (pow(IP.requil/ZZ,12) - 2*pow(IP.requil/ZZ,6));
}

double WL (double zz, programInputs &IP)
{
  return 1.0/(1.0+exp(IP.alpha*(zz-IP.zl)));
}

double WR (double zz, programInputs &IP)
{
  return 1.0/(1.0+exp(-1.0*IP.alpha*(zz-IP.zr)));
}

double VVac (double zz, programInputs &IP)
{
  double dist = IP.zr - IP.zl;
  return IP.echarge*IP.bias*(1.0 - (zz - IP.zl)/dist)/27.21; //output in atomic units
}

double Ve (double zz, programInputs &IP)
{
  double VR = IP.workfxn;
  double VL = IP.workfxn + (IP.echarge * IP.bias)/27.21; //in a.u.
  return WR(zz,IP) * VR + WL(zz,IP) * VL + (1.0 - WR(zz,IP) - WL(zz,IP))*VVac(zz,IP);
}

double VeMS (double ZZ, programInputs &IP)
{
  if (IP.eaffin == 0)
    return 0.0;
  else
    return IP.eaffin - (IP.cconst/(27.21*(ZZ-IP.zl)));
}

double Wcoup (double zz, double ZZ, programInputs &IP)
{
  if (IP.cconst == 0)
    return exp(-1.0*(zz-ZZ)*(zz-ZZ)/(2.0*pow(IP.sigma,2)));
  else
  {
    double diff   = abs(zz-ZZ);
    double parab  = (diff/IP.rad) * (diff/IP.rad) - 1.0;
    double on     = 1.0/(1.0 + exp(IP.del*(diff-IP.rad)));
    double off    = 1.0/(1.0 + exp(-1.0*IP.del*(diff+IP.rad)));
    return -1.0*parab*on*off;
  }
}

double ElecParab (double zz, double ZZ, programInputs &IP)
{
  double diff = abs(zz-ZZ);
  return pow((diff/IP.rad),2)*abs(IP.eaffin);
}

double VeM (double zz, double ZZ, programInputs &IP)
{
  return Wcoup(zz,ZZ,IP)*(VeMS(ZZ,IP) - Ve(zz,IP));
}

double V0 (double zz, double ZZ, programInputs &IP)
{
  return VM(ZZ,IP) + Ve(zz,IP);
}

double absorbingPotential(double XX, programInputs &IP)
{
  double c = 2.622;
  double xx;
  if (XX <= (IP.xmin + IP.lbarrier))
  {
    xx = -1.0 * (c +0.001)*(XX-(IP.xmin + IP.lbarrier))/IP.lbarrier;
    return IP.lbarrieritn*(4.0/((c-xx)*(c-xx))+4.0/((c+xx)*(c+xx))-8.0/(c*c))*(cosh(log(2.0+sqrt(3.0))*xx/c)-1.0);
  }
  else if (XX >= (IP.xmax - IP.rbarrier))
  {
    xx = (c+0.001)*(XX-(IP.xmax - IP.rbarrier))/IP.rbarrier;
    return IP.rbarrieritn*(4.0/((c-xx)*(c-xx))+4.0/((c+xx)*(c+xx))-8.0/(c*c))*(cosh(log(2.0+sqrt(3.0))*xx/c)-1.0);
  }
  else
    return 0;
}

//  Absorbing potential in the molecular dimension
double absorbingPotentialMol(double XX, programInputs &IP)
{
  double c = 2.622;
  double xx;
  double Lbarrier = 3.0;
  double Rbarrier = 5.0;
  if (XX <= (IP.ymin + Lbarrier))
  {
    xx = -1.0 * (c +0.001)*(XX-(IP.ymin + Lbarrier))/Lbarrier;
    return (2*M_PI*M_PI/(Lbarrier*Lbarrier))*(4.0/((c-xx)*(c-xx))+4.0/((c+xx)*(c+xx))-8.0/(c*c))*(cosh(log(2.0+sqrt(3.0))*xx/c)-1.0);
  }
  else if (XX >= (IP.ymax - Rbarrier))
  {
    xx = (c+0.001)*(XX-(IP.ymax - Rbarrier))/Rbarrier;
    return (2*M_PI*M_PI/(Rbarrier*Rbarrier))*(4.0/((c-xx)*(c-xx))+4.0/((c+xx)*(c+xx))-8.0/(c*c))*(cosh(log(2.0+sqrt(3.0))*xx/c)-1.0);
  }
  else
    return 0;
}

double HPotentialNoAbs (double zz, double ZZ, programInputs &IP)
{
  return (VM(ZZ,IP) + Ve(zz,IP) + VeM(zz,ZZ,IP));
}

cplx HPotential (double zz, double ZZ, programInputs &IP)
{
  return cplx(0.0,-1.0)*(absorbingPotential(zz,IP) + absorbingPotentialMol(ZZ,IP))+ HPotentialNoAbs(zz,ZZ,IP);
}

double wvfxnNuclear (double ZZ, programInputs &IP)
{
  return hoWvfxn(0, ZZ, IP.requil,IP.omega,IP.m_C60,IP.hbar);
}

cplx wvfxnElectron (double zz, programInputs &IP)
{
  double coeff = sqrt(sqrt(1.0/(M_PI*IP.xi*IP.xi)));
  cplx arg1    = cplx(0.0,1.0) * IP.elecvel * zz / IP.hbar;
  arg1         = exp(arg1);
  double arg2  = exp(-1.0 * (zz-IP.elecpos) * (zz-IP.elecpos)/(2.0*IP.xi*IP.xi));
  cplx val     = coeff * arg1 * arg2;
  return val;
}

double hermite(int n, double x)
{
  if (n == 0)
    return 1;
  else if (n == 1)
    return 2*x;
  else
    return 2*x*hermite(n-1,x) - 2*(n-1)*hermite(n-2,x);
}

double factorial(int n)
{
  double fact;
  if (n==0)
    fact = 1;
  else
    fact = n*factorial(n-1);
  return fact;
}

double hoWvfxn(int n, double x, double xcen, double omega, double m, double hbar)
{ //nth harmonic oscillator solution evaluated at x, centered at xcen
  double coeff1 = 1.0/sqrt(pow(2,n)*factorial(n));
  double coeff2 = sqrt(sqrt(m*omega/M_PI/hbar));
  double coeff  = coeff1 * coeff2;
  x             -= xcen;
  double fxn    = coeff * exp(-1*m*omega*x*x/2/hbar) * hermite(n,sqrt(m*omega/hbar)*x);
  return fxn;
}

int kron_delta(int x, int y)
{
  if (x == y)
    return 1;
  else
    return 0;
}
