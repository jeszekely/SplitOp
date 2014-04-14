#include <memory>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <random>
#include <stdexcept>

#include "vector.hpp"
#include "matrix.hpp"
#include "utilities.hpp"

using namespace std;

vectorMatrix::vectorMatrix(const int nr, const int nc) : matrixReal(nr,nc){}
vectorMatrix::vectorMatrix(const vectorMatrix& o) : matrixReal(o){}
vectorMatrix::vectorMatrix(vectorMatrix&& o) : matrixReal(move(o)){}
vectorMatrix::vectorMatrix(const matrixReal& o) : matrixReal(o){}
vectorMatrix::vectorMatrix(matrixReal&& o) : matrixReal(move(o)){}

void vectorMatrix::normalize(const int n, const double factor)
{
  double norm = sqrt(dot(n,n));
  scaleVec(n,factor/norm);
}

void vectorMatrix::normalize(const int n)
{
  normalize(n,1.0);
}

void vectorMatrix::normalizeAll()
{
  for (int ii = 0; ii<ncols; ii++)
    normalize(ii);
}

//Gram-Schmidt orthonormalization procedure
void vectorMatrix::orthonorm(const int n)
{
  for (int ii = 0; ii<ncols; ii++)
  {
    if (ii == n) continue; //skips orthogonalization with self
    double factor = dot(ii,n);
    setSub(0,n,*getSub(0,n,nrows,1)-*getSub(0,ii,nrows,1)*factor);
  }
  normalize(n);
}

void vectorMatrix::orthonormAll()
{
  for (int ii = 1; ii < ncols; ii++)
   {
    for (int jj = 0; jj < ii; jj++)
     {
        double factor = dot(ii,jj);
        setSub(0,ii,*getSub(0,ii,nrows,1) -*getSub(0,jj,nrows,1)*factor);
     }
     double norm =  dot(ii,ii);
     scaleVec(ii,1.0/sqrt(norm));
   }
   normalizeAll();
 }

shared_ptr<vectorMatrix> vectorMatrix::canonical_orthogonalization(const double thresh) const {
  matrixReal S(*this | *this);
  vector<double> S_eigs(nc(), 0.0);
  S.diagonalize(S_eigs.data());

  const int out_dim = count_if(S_eigs.begin(), S_eigs.end(), [&thresh] (const double v) { return (v > thresh); });
  auto tmp = make_shared<vectorMatrix>(S.nr(), out_dim);
  for (int i = 0, current = 0; i < nc(); ++i) {
    if (S_eigs[i] > thresh) {
      daxpy_(S.nr(), 1.0/std::sqrt(S_eigs[i]), &S(0,i), 1, &tmp->element(0,current++), 1);
    }
  }

  return make_shared<vectorMatrix>(*this * *tmp);
}

//vector dot product where vectors are two columns in matrix
double vectorMatrix::dot(const int ii, const int jj) const
{
  assert(ncols > max(ii,jj));
  return ddot_(nrows, &element(0,ii), 1, &element(0,jj), 1);
}

//vector dot product where vectors are two columns in separate matricies
double vectorMatrix::dot(const vectorMatrix& o, const int ii, const int jj) const
{
  assert(ncols > ii);
  assert(o.ncols > jj);
  assert(nrows == o.nrows);
  return ddot_(nrows, &element(0,ii), 1, &o.element(0,jj), 1);
}

matrixReal vectorMatrix::vec(const int vecN)
{
  matrixReal out(*getSub(0,vecN,nrows,1));
  return out;
}

void vectorMatrix::diagonalizeSub(double* eigVals,int nr, int nc)
{
  assert (nrows == ncols);
  int info;
  int lwork = -1;
  double wkopt;
  dsyev_("V", "U", nr, data(), nrows, eigVals, &wkopt, lwork, info);
  lwork = int(wkopt);
  std::unique_ptr <double[]> work (new double [lwork]);
  dsyev_("V", "U", nr, data(), nrows, eigVals, work.get(), lwork, info);
  if (info > 0)
    throw std::runtime_error("Unable to diagonalize matrix");
  return;
}

