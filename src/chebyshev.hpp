#ifndef SPLITOP_CHEBYSHEV
#define SPLITOP_CHEBYSHEV

#include "wvfxn.hpp"
#include "input_parser.hpp"
#include <boost/math/special_functions/bessel.hpp>
#include "fftw3.h"
#include <cmath>
#include <complex>

class Chebyshev1D
{
protected:
  size_t nx;
  int nthreads;
  int polyterms;
  double Vmin,Vmax,Tmin,Tmax;

public:
  double xmin, xmax, xstep, pstep;
  double itpdt, dt, simtime, runtime;
  std::shared_ptr<Array1D<double>> xgrid, pgrid, Tgrid,Vgrid;
  std::shared_ptr<wvfxn1D> wvfxn;

  Chebyshev1D(programInputs &IP);
  void initializeTDSE(std::function<double(double)>,std::function<double(double)>);
  void propagateStep();
  void ApplyHNorm(std::shared_ptr<wvfxn1D> W);
};

double ChebyshevCoeff(int nn, double alpha);

#endif