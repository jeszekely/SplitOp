#include "wvfxn.hpp"

using namespace std;

wvfxn1D::wvfxn1D(const int Nx, const double xi, const double xs) : Array1D<cplx>(Nx,cplx(xi),cplx(xs)){mass = hbar = 1.0;}
wvfxn1D::wvfxn1D(const int Nx, const cplx xi, const cplx xs) : Array1D<cplx>(Nx,xi,xs){mass = hbar = 1.0;}
wvfxn1D::wvfxn1D(const int Nx, const double xi, const double xs, const double hb, const double m) : Array1D<cplx>(Nx,cplx(xi),cplx(xs)), hbar(hb), mass(m){}
wvfxn1D::wvfxn1D(const int Nx, const cplx xi, const cplx xs, const double hb, const double m) : Array1D<cplx>(Nx,xi,xs), hbar(hb), mass(m){}

wvfxn1D::wvfxn1D(const wvfxn1D& o) : Array1D<cplx>(o){}
wvfxn1D::wvfxn1D(wvfxn1D&& o) : Array1D<cplx>(move(o)){}

wvfxn1D& wvfxn1D::operator=(const wvfxn1D& o)
{
  assert(nx == o.nx);
  copy_n(o.data(), o.size(), data());
  return *this;
}

double wvfxn1D::getNorm()
{
  auto S = *this | *this;
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
  cplx psi_deriv = cplx(0,-1.0)*deriv_16(xx);
  return real(psi_deriv*psi_star*hbar/mass);
}

double wvfxn1D::m() {return mass;}
double wvfxn1D::hb() {return hbar;}

wvfxn2D::wvfxn2D(const int Nx, const int Ny, const double xs, const double ys) : Array2D<cplx>(Nx,Ny,cplx(xs),cplx(ys)) {mass1 = mass2 = hbar = 1.0;}
wvfxn2D::wvfxn2D(const int Nx, const int Ny, const cplx xs, const cplx ys) : Array2D<cplx>(Nx,Ny,xs,ys){mass1 = mass2 = hbar = 1.0;}
wvfxn2D::wvfxn2D(const int Nx, const int Ny, const double xs, const double ys, const double hb, const double m1, const double m2) : Array2D<cplx>(Nx,Ny,cplx(xs),cplx(ys)), hbar(hb), mass1(m1), mass2(m2){}
wvfxn2D::wvfxn2D(const int Nx, const int Ny, const cplx xs, const cplx ys, const double hb, const double m1, const double m2) : Array2D<cplx>(Nx,Ny,xs,ys), hbar(hb), mass1(m1), mass2(m2) {}
wvfxn2D::wvfxn2D(const wvfxn2D& o) : Array2D<cplx>(o){}
wvfxn2D::wvfxn2D(wvfxn2D&& o) : Array2D<cplx>(move(o)){}

wvfxn2D& wvfxn2D::operator=(const wvfxn2D& o)
{
  assert(nx == o.nx && ny == o.ny);
  copy_n(o.data(), o.size(), data());
  return *this;
}

double wvfxn2D::getNorm()
{
  auto S = *this | *this;
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
  Array1D <double> probcurr(ny,0.0,real(ystep));
  for (int ii = 0; ii < ny; ii++)
  {
    psi_star = conj(element(xx,ii));
    psi_deriv = cplx(0.0,-1.0)*deriv_16_x(xx,ii);
    probcurr(ii) = real(psi_deriv*psi_star);
  }
  return (mass1/hbar)*probcurr.integrate_rect();
}
double wvfxn2D::flux_y(const int yy)
{
  cplx psi_star, psi_deriv;
  Array1D <double> probcurr(nx,0.0,real(xstep));
  for (int ii = 0; ii < nx; ii++)
  {
    psi_star = conj(element(ii,yy));
    psi_deriv = cplx(0.0,-1.0)*deriv_16_y(ii,yy);
    probcurr(ii) = real(psi_deriv*psi_star);
  }
  return (mass2/hbar)*probcurr.integrate_rect();
}
double wvfxn2D::m1() {return mass1;}
double wvfxn2D::m2() {return mass2;}
double wvfxn2D::hb() {return hbar;}