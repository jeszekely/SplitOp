#ifndef SPLIPOP_SPLITOP
#define SPLITOP_SPLITOP

#include "wvfxn.hpp"
#include "input_parser.hpp"
#include "fftw3.h"
#include "fftw3_mkl.h"
#include <cmath>

//typedef std::complex<double> cplx;

class SplitOp1D
{
protected:
  size_t nx;
  int nthreads;
  fftw_plan forplan, backplan;
  double xmin, xmax, xstep, pstep;

public:
  double dt, simtime, runtime;
  std::shared_ptr<Array1D<double>> xgrid, pgrid, Vgrid, Tgrid;
  std::shared_ptr<Array1D<cplx>> KinetOp, PotenOp;
  std::shared_ptr<wvfxn1D> wvfxn;

  SplitOp1D(programInputs &IP);
  void initializeTDSE(std::function<cplx(double)>,std::function<double(double)>);
  // void propagateStep();
  // void propagateNSteps(int);

  // //ITP Specific functions
  // //number of states, array to place them, convergence criterion
  // void initializeITP(Array1D <double> &, Array1D <double> &);
  // void propagateITP(int,Array2D<cplx> &, double);
  // void initializeGuess();

  // //Void output functions
  // void outputWvfxn();
  // void outputOperator(Array1D &);
};

#endif
