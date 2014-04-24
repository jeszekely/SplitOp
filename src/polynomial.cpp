#include "polynomial.hpp"
#include <boost/math/special_functions/bessel.hpp>

using namespace std;

double ChebyshevCoeff(int nn, double alpha)
{
  return 2.0*boost::math::cyl_bessel_j(nn,alpha);
}

shared_ptr<polynomial<cplx>> ClenshawChebyshevProp(int nn, double alpha)
{
  std::shared_ptr<polynomial<T>> p0,p1,pn;
  polynomial<T> a1(2);
  a1(1)          = 2;
  p0             = std::make_shared<polynomial<T>>(1);
  p0->element(0) = T(1.0);
  p1             = std::make_shared<polynomial<T>>(2);
  p1->element(1) = T(1.0);
  if (nn == 0) return p0;
  if (nn == 1) return p1;
  int kk = 1;
  while (kk < nn)
  {
    kk++;
    pn = std::make_shared<polynomial<T>>(*p1*a1-(*p0));
    p0 = p1;
    p1 = pn;
  }
  return pn;
}
