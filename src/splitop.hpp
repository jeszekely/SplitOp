#ifndef SPLIPOP_SPLITOP
#define SPLITOP_SPLITOP

#include "wvfxn.hpp"
#include <fftw3.h>

//typedef std::complex<double> cplx;

class SplitOp1D
{
protected:
  size_t nx;
  double xinit, xstep, pstep;
  double dt;
  int nthreads;
  fftw_plan forplan, backplan;
public:
  Array1D <double> xgrid, pgrid;
  Array1D <cplx> KinetOp, PotenOp;
  wvfxn1D wvfxn;
  void initializeTDSE(Array1D <double> &, Array1D <double> &);
  void propagateStep();
  void propagateNSteps(int);

  //ITP Specific functions
  //number of states, array to place them, convergence criterion
  void propagateITP(int,Array2D<cplx> &, double);
  void initializeITP(Array1D <double> &, Array1D <double> &);
  void initializeGuess();
}

#endif
