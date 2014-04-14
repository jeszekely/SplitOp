#ifndef DMRG_VECTOR
#define DMRG_VECTOR

#include <memory>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <random>

#include "matrix.hpp"
#include "utilities.hpp"

class vectorMatrix : public matrixReal
{
public:
  vectorMatrix(const int nr, const int nc);
  vectorMatrix(const vectorMatrix&);
  vectorMatrix(vectorMatrix&&);
  vectorMatrix(const matrixReal&);
  vectorMatrix(matrixReal&&);

//Vector-Vector Operations
  double dot(const int ii, const int jj) const;
  double dot(const vectorMatrix& o, const int ii, const int jj) const;

//  Apply a scaling factor to one vector
  void scaleVec(const int n, const double factor)
  {
    std::for_each(&element(0,n), &element(0,n)+nrows, [&factor](double & p){p*=factor;});
  }

  void normalize(const int n, const double scale);
  void normalize(const int n);
  void normalizeAll();

  void orthonorm(int n);
  void orthonormAll();

  std::shared_ptr<vectorMatrix> canonical_orthogonalization(const double thresh = 1.0e-8) const;

  matrixReal vec(const int vecN);
  void diagonalizeSub(double* eigVals,int nr, int nc);

};
#endif
