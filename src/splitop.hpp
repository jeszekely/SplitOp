#ifndef SPLIPOP_SPLITOP
#define SPLITOP_SPLITOP

#include "wvfxn.hpp"
#include "input_parser.hpp"
#include "fftw3.h"
#include <cmath>
#include <complex>

class SplitOp1D
{
protected:
  size_t nx;
  int nthreads;
  fftw_plan forplan, backplan;
  double xmin, xmax, xstep, pstep;

public:
  double itpdt, dt, simtime, runtime;
  std::shared_ptr<Array1D<double>> xgrid, pgrid, Tgrid;
  std::shared_ptr<Array1D<cplx>> Vgrid, KinetOp, PotenOp;
  std::shared_ptr<wvfxn1D> wvfxn;

  SplitOp1D(programInputs &IP);
  ~SplitOp1D();
  void initializeTDSE(std::function<cplx(double)>,std::function<double(double)>);
  void propagateStep();
  void propagateNSteps(int);

  // //ITP Specific functions
  // //number of states, array to place them, convergence criterion
  void initializeITP(std::function<cplx(double)>,std::function<double(double)>);
  void propagateITP(Array2D<cplx> &, double);
  void updateITP();
  void initializeGuess(std::vector<std::shared_ptr<wvfxn1D>>);

  // //Void output functions
  // void outputWvfxn();
  // void outputOperator(Array1D &);
};

// class SplitOp2D
// {
// protected:
//   size_t nx, ny;
//   int nthreads;
//   fftw_plan forplan, backplan;
//   double xmin, xmax, xstep, pstep;
//   double ymin, ymax, ystep, qstep;

// public:
//   double dt, simtime, runtime;
//   std::shared_ptr<Array1D<double>> xgrid, pgrid, ygrid, qgrid;
//   std::shared_ptr<Array2D<double>> Tgrid;
//   std::shared_ptr<Array2D<cplx>> Vgrid, KinetOp, PotenOp;
//   std::shared_ptr<wvfxn2D> wvfxn;

//   SplitOp1D(programInputs &IP);
//   ~SplitOp1D();
//   void initializeTDSE(std::function<cplx(double)>,std::function<double(double)>);
//   void propagateStep();
//   void propagateNSteps(int);

//   // //Void output functions
//   // void outputWvfxn();
//   // void outputOperator(Array1D &);
// };

#endif
