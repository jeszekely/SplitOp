#include "polynomial.hpp"
#include <boost/math/special_functions/bessel.hpp>

using namespace std;

double ChebyshevCoeff(int nn, double alpha)
{
  if (nn == 0) return boost::math::cyl_bessel_j(nn,alpha);
  else return 2.0*boost::math::cyl_bessel_j(nn,alpha);
}

shared_ptr<polynomial<cplx>> ClenshawChebyshevProp(int nn, double a)
{
  assert(nn >= 2);
  polynomial<cplx> phi0(1);
  phi0(0)  = 1;
  polynomial<cplx> phi1(2);
  phi1(1)  = cplx(0.0,-1.0);
  polynomial<cplx> bkp1(1);
  polynomial<cplx> bkp2(1);
  polynomial<cplx> alpha(2);
  alpha(1) = cplx(0.0,-2.0);
  polynomial<cplx> beta(1);
  beta(0)  = 1; //-1 for normal chebyshev iterative definition
  int kk   = nn;
  polynomial<cplx> bk(1);
  polynomial<cplx> ak(1);
  polynomial<cplx> a0(1);
  while (true)
  {
    ak(0) = ChebyshevCoeff(kk,a);
    bk    = alpha*bkp1 + beta*bkp2 + ak;
    //cout << kk << ":" << endl << bk << bkp1 <<  bkp2 << phi1*bkp1+beta*phi0*bkp2 << endl << endl;
    if (kk == 0) break;
    bkp2  = bkp1;
    bkp1  = bk;
    --kk;
  }
  shared_ptr<polynomial<cplx>> S = make_shared<polynomial<cplx>>(phi1*bkp1+beta*phi0*bkp2+phi0*ak);
  return S;
}


