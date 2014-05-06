#ifndef DMRG_MATRIX
#define DMRG_MATRIX

#include <algorithm>
#include <assert.h>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <complex>
#include "utilities.hpp"

typedef std::complex<double> cplx;

template <typename T> class matrixBase
{
protected:
    static unsigned int memSize;
    size_t nrows, ncols;
    std::unique_ptr<T[]> vals;

public:
//  Constructor
  matrixBase(const int nr, const int nc) : nrows(nr), ncols(nc), vals(std::unique_ptr<T[]>(new T[nc*nr]))
  {
    zero();
    memSize += sizeof(T)*size(); //number of bytes allocated to instance of class
  }

//  Copy Constructor
  matrixBase(const matrixBase& o) : nrows(o.nrows), ncols(o.ncols), vals(std::unique_ptr<T[]>(new T[nrows*ncols]))
  {
    std::copy_n(o.vals.get(), nrows*ncols, vals.get());
    memSize += sizeof(T)*size(); //number of bytes allocated to instance of class
  }

//  Move Constructor
  matrixBase(matrixBase&& o) : nrows(o.nrows), ncols(o.ncols), vals(std::move(o.vals)) { o.nrows = 0; o.ncols = 0; };

//  Destructor
  ~matrixBase(){memSize -= sizeof(T)*size();} //number of bytes allocated to instance of class

//  Access functions
  size_t size() const { return nrows * ncols; }

  T* data() { return vals.get(); }
  const T* data() const { return vals.get(); }

//  Fill with zeroes
  void zero()
  {
    std::fill_n(vals.get(), nrows*ncols, T(0.0));
  }

  size_t nr() const {return nrows;}
  size_t nc() const {return ncols;}

  void random()
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> dis(-2, 2);
    std::generate_n(vals.get(), nrows*ncols, [&dis, &gen](){return dis(gen);});
  }

//  Accessor functions
  T& element(const int row, const int col)
  {
    return vals[col*nrows + row];
  }

  const T& element(const int row, const int col) const
  {
    return vals[col*nrows + row];
  }

  T& operator()(const int row, const int col)
  {
    return element(row,col);
  }

  const T& operator()(const int row, const int col) const
  {
    return vals[col*nrows + row];
  }

//  Set only diagonal elements to unity
  void makeIdentity()
  {
    zero();
    for (int ii = 0; ii < std::min(ncols,nrows); ii++)
      vals[ii+ii*nrows] = T(1.0);
  }

//  Compute the trace
  T trace()
  {
    T sum = T(0.0);
    for (int ii = 0; ii < std::min(int(ncols),int(nrows)); ii++)
      sum += element(ii,ii);
    return sum;
  }

//  Apply a scaling factor to all elements
  void scale(const T a)
  {
    std::for_each(data(), data()+size(), [&a](T& p){p*=a;});
  }

//  Print memory usage for all matrices
  void printMem() const
  {
    std::cout << "Current memory allocated to this matrix: " << size()*sizeof(T) << " bytes." << std::endl;
    std::cout << "Total memory allocated for matrix storage: " << memSize << " bytes." << std::endl;
    return;
  }

//Extract a portion of the matrix starting at element(r,c), get (nr x nc matrix)
  template <class U>
  std::shared_ptr<U> getSub_impl(int r, int c, int nr, int nc) const
  {
    assert(r + nr <= nrows && c + nc <= ncols);
    auto out = std::make_shared<U>(nr,nc);
    for (int jj = 0; jj < nc; jj++)
      std::copy_n(&element(r,c+jj),nr,&out->element(0,jj));
    return out;
 }

//Place matrix o at position (r,c)
  template <typename U>
  void setSub(int r, int c, U o)
  {
    assert(r + o.nrows <= nrows && c + o.ncols <= ncols);
    for (int jj = 0; jj < o.ncols; jj++)
      std::copy_n(&o(0,jj), o.nrows, &element(r,c+jj));
    return;
  }

  template <typename U> friend std::ostream &operator<<(std::ostream &out, const matrixBase <U> &o);
};

class matrixReal : public matrixBase<double>
{
public:
  matrixReal(const int nr, const int nc);
  matrixReal(const matrixReal&);
  matrixReal(matrixReal&&);

//Matrix-Matrix operations
  matrixReal& operator=(const matrixReal&);
  matrixReal operator*(const matrixReal&) const;
  matrixReal& operator*=(const matrixReal&);
  matrixReal operator+(const matrixReal&) const;
  matrixReal& operator+=(const matrixReal&);
  matrixReal operator-(const matrixReal&) const;
  matrixReal& operator-=(const matrixReal&);
  matrixReal operator|(const matrixReal&) const;
  matrixReal operator^(const matrixReal&) const;

//Scalar-Matrix Operations
//Note: binary scalar operations only work as rhs operators at the moment
  matrixReal operator*(const double&) const;
  matrixReal operator/(const double&) const;
  matrixReal& operator*=(const double&);
  matrixReal& operator/=(const double&);

  double dot_product(const matrixReal& o) const;
  double norm() const;
  double rms() const;
  double variance() const;

//Vector-Vector Operations
  double operator%(const matrixReal& o) const;

//BLAS and LAPACK routines
//Diagonalize matrix, place eigenvalues in the vector prodived
//NOTE: Assumes a symmetric matrix
  void diagonalize(double* eigVals);
  void diagonalize(double* eigVals, bool getLowEigVal, int keepNum);
  std::shared_ptr<matrixReal> transpose() const;
  std::tuple<std::shared_ptr<matrixReal>, std::shared_ptr<matrixReal>>svd(std::vector<double>&);

//Compute the kronecker product of two matrices
//Note: This will be used as a shortcut to create a superblock matrix for small site matrices
//This should NOT be used for larger systems as it is very memory intensive
  matrixReal kron(matrixReal &o) const;

  std::shared_ptr<matrixReal> getSub(int ii, int jj, int kk, int ll) const
  {
    return getSub_impl<matrixReal>(ii,jj,kk,ll);
  }

  void ax_plus_y(const double a, matrixReal &o);
};

//Overload the << operator to print a matrix
template <typename T>
std::ostream &operator<<(std::ostream &out, const matrixBase <T> &o)
{
  for (int row = 0; row < std::min(10,int(o.nrows)); row++)
  {
    for (int col = 0; col < std::min(10,int(o.ncols)); col++)
    {
      out << std::setprecision(3) << o(row,col) << "\t";
    }
    out << "\n";
  }
  out << std::endl;
  return out;
};

template <typename T> unsigned int matrixBase<T>::memSize = 0;



//class matrixComplex

class matrixComp : public matrixBase<cplx>
{
public:
  matrixComp(const int nr, const int nc);
  matrixComp(const matrixComp&);
  matrixComp(matrixComp&&);

//Matrix-Matrix operations
  matrixComp& operator=(const matrixComp&);
  matrixComp operator*(const matrixComp&) const;
  matrixComp& operator*=(const matrixComp&);
  matrixComp operator+(const matrixComp&) const;
  matrixComp& operator+=(const matrixComp&);
  matrixComp operator-(const matrixComp&) const;
  matrixComp& operator-=(const matrixComp&);
  // matrixComp operator|(const matrixComp&) const;
  // matrixComp operator^(const matrixComp&) const;
  void getEigvals(double* eigVals);


//Scalar-Matrix Operations
//Note: binary scalar operations only work as rhs operators at the moment
  matrixComp operator*(const cplx&) const;
  matrixComp operator/(const cplx&) const;
  matrixComp& operator*=(const cplx&);
  matrixComp& operator/=(const cplx&);

  // double dot_product(const matrixComp& o) const;
  // double norm() const;
  // double rms() const;
  // double variance() const;

};


#endif
