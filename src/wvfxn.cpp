#include "wvfxn.hpp"

typedef std::complex<double> cplx;

using namespace std;

wvfxn1D::wvfxn1D(const int Nx, const double xi, const double xs) : Array1D<cplx>(Nx,xi,xs){mass = hbar = 1.0;}
wvfxn1D::wvfxn1D(const wvfxn1D& o) : Array1D<cplx>(o){}
wvfxn1D::wvfxn1D(wvfxn1D&& o) : Array1D<cplx>(move(o)){}

double wvfxn1D::getNorm()
{
  double S = *this | *this;
  return real(S.integrate_rect());
}
void wvfxn1D::normalize()
{
  double norm = getNorm();
  assert (norm > 0);
  scale(1.0/sqrt(norm));
}

double wvfxn1D::flux(const int xx)
{
  cplx psi_star = conj(vals[xx]);
  cplx psi_deriv = -1.0j*deriv_16(xx);
  return real(psi_deriv*psi_star*hbar/mass);
}
/*
wvfxn2D::wvfxn2D(const int Nx, const int Ny, const double xs, const double ys) : Array1D<cplx>(Nx,Ny,xs,ys){mass = hbar = 1.0;}
wvfxn2D::wvfxn2D(const wvfxn2D& o) : Array2D<cplx>(o){}
wvfxn2D::wvfxn2D(wvfxn2D&& o) : Array2D<cplx>(move(o)){}

double wvfxn2D::getNorm()
{
  double S = *this | *this;
  return real(S.integrate_rect());
}

void wvfxn2D::normalize()
{
  double norm = getNorm();
  assert (norm > 0);
  scale(1.0/sqrt(norm));
}

double wvfxn2D::flux_x(const int xx)
{
  cplx psi_star, psi_deriv;
  Array1D <double> probcurr(ny,yinit,ystep);
  for (int ii = 0; ii < ny; ii++)
  {
    psi_star = conj(element(xx,ii));
    psi_deriv = -1.0j*deriv_16_x(xx,ii);
    probcurr(ii) = real(psi_deriv*psi_star);
  }
  return (mass/hbar)*probcurr.integrate_rect();
}
double wvfxn2D::flux_y(const int yy)
{
  cplx psi_star, psi_deriv;
  Array1D <double> probcurr(nx,xinit,xstep);
  for (int ii = 0; ii < nx; ii++)
  {
    psi_star = conj(element(ii,yy));
    psi_deriv = -1.0j*deriv_16_y(ii,yy);
    probcurr(ii) = real(psi_deriv*psi_star);
  }
  return (mass/hbar)*probcurr.integrate_rect();
}

*/