#include "polynomial.hpp"
#include <boost/math/special_functions/bessel.hpp>

using namespace std;

double ChebyshevCoeff(int nn, double alpha)
{
  return 2.0*boost::math::cyl_bessel_j(nn,alpha);
}

shared_ptr<polynomial<cplx>> ClenshawChebyshevProp(int nn, double a)
{
  assert(nn > 2);
  polynomial<cplx> phi1(1);
  phi1(0)  = 1;
  polynomial<cplx> phi2(2);
  phi2(1)  = 1;
  polynomial<cplx> bnp1(1);
  polynomial<cplx> bnp2(1);
  polynomial<cplx> alpha(2);
  alpha(1) = 2;
  polynomial<cplx> beta(1);
  beta(0)  = -1;
  int kk   = nn;
  polynomial<cplx> bn(1);

  while (kk > 1)
  {
    // cout << kk << ":" << endl << bnp2 << bnp1 << bn << endl;
    bnp2  = bnp1;
    bnp1  = bn;
    cout << bn;
    polynomial<cplx> ak(1);
    ak(0) = 1.0;//pow(cplx(0.0,1.0),kk)* ChebyshevCoeff(kk,a);
    bn    = alpha*bnp1 + beta*bnp2 + ak;
    kk--;
  }
  shared_ptr<polynomial<cplx>> S = make_shared<polynomial<cplx>>(phi2*bnp1+beta*phi1*bnp2);
  return S;
}


