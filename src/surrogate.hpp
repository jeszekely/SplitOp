#ifndef SPLITOP_SURROGATE
#define SPLITOP_SURROGATE

#include <complex>
#include <assert.h>
#include "arrays.hpp"
#include "input_parser.hpp"
#include "wvfxn.hpp"

typedef std::complex<double> cplx;

//Wavefuction array class for all 1D functions needed
class wvfxn1DArray : public Array2D<cplx>
{
public:
  double hbar;
  double mass;

  //number of wavefunctions, input parameters (hbar, mass)
  wvfxn1DArray(const int nn, const double xstep, programInputs &IP);

  //wvfxn1DArray(const wvfxn1DArray&);
  //wvfxn1DArray(wvfxn1DArray&&);
  //wvfxn1DArray& operator=(const wvfxn1DArray&);

//Wavefunction specific operations
  double getNorm(const int nn);
  void normalize(const int nn);
  double flux(const int nn, const int xx);
  double hb();
  double m();
};



class SplitOp1DArray
{
protected:
  size_t nx;
  int nthreads;
  std::vector<fftw_plan> forplanArray, backplanArray;
  std::vector<std::make_shared<Array2D<cplx>>> HInteractions;

public:
  double xmin, xmax, xstep, pstep;
  double dt, simtime, runtime;
  std::shared_ptr<Array1D<double>> xgrid, pgrid, Tgrid;
  std::shared_ptr<Array1D<cplx>> Vgrid, KinetOp, PotenOp;
  std::shared_ptr<wvfxn1DArray> wvfxn;

  SplitOp1DArray(programInputs &IP);
  ~SplitOp1DArray();

  void propagateStep(int ii = 1);

};

#endif