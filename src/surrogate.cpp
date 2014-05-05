#include "surrogate.hpp"

using namespace std;

wvfxn1DArray::wvfxn1DArray(const int nn, const double xstep, programInputs &IP) : Array2D<cplx>(IP.nx,nn,xstep,1.0)
{
  mass = IP.m_electron;
  hbar = IP.hbar;
}

double wvfxn1DArray::getNorm(const int nn)
{
  wvfxn1D S(nx,0.0,xstep);
  copy_n(&element(0,nn),nx,S.data());
  auto T = S | S;
  return real(T.integrate_rect());
}

void wvfxn1DArray::normalize(const int nn)
{
  double norm = getNorm(nn);
  assert (norm > 0);
  double fact = 1.0/sqrt(norm);
  std::for_each(&element(0,nn),&element(0,nn)+nx, [&fact](cplx &p){p*=fact;});
}

double wvfxn1DArray::flux(const int nn, const int xx)
{
  cplx psi_star = conj(element(nn,xx));
  cplx psi_deriv = cplx(0,-1.0)*deriv_16_x(nn,xx);
  return real(psi_deriv*psi_star*hbar/mass);
}

double wvfxn1DArray::m() {return mass;}
double wvfxn1DArray::hb() {return hbar;}