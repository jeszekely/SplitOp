#include "matrix.hpp"

using namespace std;

matrixReal::matrixReal(const int nr, const int nc) : matrixBase<double>(nr,nc){}
matrixReal::matrixReal(const matrixReal& o) : matrixBase<double>(o){}
matrixReal::matrixReal(matrixReal&& o ) : matrixBase<double>(std::move(o)){}

matrixReal& matrixReal::operator=(const matrixReal& o)
{
  assert(nrows == o.nrows && ncols == o.ncols);
  copy_n(o.data(), o.size(), data());
  return *this;
}

matrixReal& matrixReal::operator*=(const matrixReal& o)
{
  *this = *this * o;
  return *this;
}

matrixReal matrixReal::operator*(const matrixReal& o) const
{
  assert(ncols == o.nrows);
  matrixReal out(nrows, o.ncols);
  dgemm_("N","N", nrows, o.ncols, o.nrows, 1.0, data(), nrows, o.data(), o.nrows, 0.0, out.data(), nrows);
  return out;
}

//vector dot product
double matrixReal::operator%(const matrixReal& o) const
{
  assert(ncols == 1 && o.ncols == 1 && nrows == o.nrows);
  return ddot_(nrows, data(), 1, o.data(), 1);
}

//Matrix multiplication, left matrix is transposed
matrixReal matrixReal::operator|(const matrixReal& o) const
{
  assert(nrows == o.nrows);
  matrixReal out(ncols, o.ncols);
  dgemm_("T","N", ncols, o.ncols, o.nrows, 1.0, data(), nrows, o.data(), o.nrows, 0.0, out.data(), ncols);
  return out;
}

//Matrix multiplication, right matrix is transposed
matrixReal matrixReal::operator^(const matrixReal& o) const
{
  assert(ncols == o.ncols);
  matrixReal out(nrows, o.nrows);
  dgemm_("N","T", nrows, o.nrows, o.nrows, 1.0, data(), nrows, o.data(), o.nrows, 0.0, out.data(), nrows);
  return out;
}

matrixReal matrixReal:: operator+(const matrixReal& o) const
{
  assert(ncols == o.ncols && nrows == o.nrows);
  matrixReal out(nrows, o.ncols);
  transform(data(), data()+size(), o.data(), out.data(), plus<double>());
  return out;
}

matrixReal& matrixReal::operator+=(const matrixReal& o)
{
  *this = *this + o;
  return *this;
}

matrixReal matrixReal::operator-(const matrixReal& o) const
{
  assert(ncols == o.ncols && nrows == o.nrows);
  matrixReal out(nrows, o.ncols);
  transform(data(), data()+size(), o.data(), out.data(), minus<double>());
  return out;
}

matrixReal& matrixReal::operator-=(const matrixReal& o)
{
  *this = *this - o;
  return *this;
}

matrixReal matrixReal::operator*(const double& a) const
{
  matrixReal out(*this);
  out *= a;
  return out;
}

matrixReal& matrixReal::operator*=(const double& a)
{
  scale(a);
  return *this;
}

matrixReal matrixReal::operator/(const double& a) const
{
  matrixReal out(*this);
  out /= a;
  return out;
}

matrixReal& matrixReal::operator/=(const double& a)
{
  scale(1.0/a);
  return *this;
}

double matrixReal::dot_product(const matrixReal& o) const {
  assert(size() == o.size());
  return ddot_(o.size(), data(), 1, o.data(), 1);
}

double matrixReal::norm() const {
  return std::sqrt(dot_product(*this));
}

double matrixReal::rms() const {
  return std::sqrt(dot_product(*this)/static_cast<double>(size()));
}

double matrixReal::variance() const {
  return dot_product(*this)/static_cast<double>(size());
}

void matrixReal::diagonalize(double* eigVals)
{
  assert (nrows == ncols);
  int info;
  int lwork = -1;
  double wkopt;
  dsyev_("V", "U", nrows, data(), nrows, eigVals, &wkopt, lwork, info);
  lwork = int(wkopt);
  std::unique_ptr <double[]> work (new double [lwork]);
  dsyev_("V", "U", nrows, data(), nrows, eigVals, work.get(), lwork, info);
  if (info > 0)
    throw std::runtime_error("Unable to diagonalize matrix");
  return;
}


void matrixReal::diagonalize(double* eigVals, bool getLowEigVal, int keepNum)
{
  assert (nrows == ncols);
  int info,iwork;
  double work;

  int lwork = -1;
  int liwork = -1;
  double abstol = 1.0e-8;
  int il, iu;
  matrixReal workEigVecs(nrows,keepNum);
  vector<int> isuppz(std::max(1,keepNum),0);
  if (getLowEigVal)
  {
      il = 1;
      iu = keepNum;
  }
  else
  {
      il = nrows - keepNum + 1;
      iu = nrows;
  }
  dsyevr_("V", "I", "U", nrows, data(), nrows, -20394857.0, 2345.2, il, iu, abstol, keepNum, eigVals, workEigVecs.data(), nrows, isuppz.data(), &work, lwork, &iwork, liwork, info);
  lwork = work;

  std::unique_ptr <double[]> workArray(new double [int(work)]);
  liwork = iwork;
  std::unique_ptr <int[]> iworkArray(new int [iwork]);

  dsyevr_("V", "I", "U", nrows, data(), nrows, -20394857.0, 2345.2, il, iu, abstol, keepNum, eigVals, workEigVecs.data(), nrows, isuppz.data(), workArray.get(), lwork, iworkArray.get(), liwork, info);
  if (info > 0)
    throw std::runtime_error("Unable to diagonalize matrix with dyevr_");
  zero();
  setSub(0,0,workEigVecs);
  return;
}

std::shared_ptr<matrixReal> matrixReal::transpose() const
{
  auto out = make_shared<matrixReal> (ncols,nrows);
  mkl_domatcopy_("C","T",nrows,ncols,1.0,data(),nrows,out->data(),ncols);
  return out;
}

tuple<shared_ptr<matrixReal>, shared_ptr<matrixReal>> matrixReal::svd(vector<double>& s)
{
  assert(s.size() >= std::min(ncols,nrows));
  auto u = make_shared<matrixReal>(nrows, nrows);
  auto vT = make_shared<matrixReal>(ncols, ncols);

  int lwork = -1;
  int info  = 0;
  dgesvd_("A","A", nrows, ncols, data(), nrows, s.data(), u->data(), nrows, vT->data(), ncols, s.data(), lwork, info);
  lwork     = s[0];
  if (lwork <= 0)
      throw runtime_error("dgesvd failed allocating lwork value");
  unique_ptr<double[]> work(new double[lwork]);
  dgesvd_("A","A", nrows, ncols, data(), nrows, s.data(), u->data(), nrows, vT->data(), ncols, work.get(), lwork, info);
  if (info != 0)
      throw runtime_error("dgesvd faied matrix decomposition");
  return make_tuple(u,vT);
}

matrixReal matrixReal::kron(matrixReal &o) const
{
  matrixReal out(nrows*o.nrows,ncols*o.ncols);
  for (int ii = 0; ii < nrows; ii++)
  {
    for (int jj = 0; jj < ncols; jj++)
    {
      out.setSub(ii*o.ncols,jj*o.nrows,o*element(ii,jj));
    }
  }
  return out;
}

void matrixReal::ax_plus_y(const double a, matrixReal &o)
{
  assert(nrows*ncols == o.nrows*o.ncols);
  daxpy_(size(), a, data(), 1, o.data(), 1);
}

  //  Hermitian Conjugate
  //  *sparsify

/*************************************
Complex Matrix Class
**************************************/

matrixComp::matrixComp(const int nr, const int nc) : matrixBase<cplx>(nr,nc){}
matrixComp::matrixComp(const matrixComp& o) : matrixBase<cplx>(o){}
matrixComp::matrixComp(matrixComp&& o ) : matrixBase<cplx>(std::move(o)){}

matrixComp& matrixComp::operator=(const matrixComp& o)
{
  assert(nrows == o.nrows && ncols == o.ncols);
  copy_n(o.data(), o.size(), data());
  return *this;
}

matrixComp& matrixComp::operator*=(const matrixComp& o)
{
  *this = *this * o;
  return *this;
}

matrixComp matrixComp::operator*(const matrixComp& o) const
{
  assert(ncols == o.nrows);
  matrixComp out(nrows, o.ncols);
  zgemm3m_("N","N", nrows, o.ncols, o.nrows, cplx(1.0), data(), nrows, o.data(), o.nrows, cplx(0.0), out.data(), nrows);
  return out;
}

// //vector dot product
// double matrixReal::operator%(const matrixReal& o) const
// {
//   assert(ncols == 1 && o.ncols == 1 && nrows == o.nrows);
//   return ddot_(nrows, data(), 1, o.data(), 1);
// }

//Matrix multiplication, left matrix is transposed
// matrixComp matrixComp::operator|(const matrixComp& o) const
// {
//   assert(nrows == o.nrows);
//   matrixComp out(ncols, o.ncols);
//   dgemm_("T","N", ncols, o.ncols, o.nrows, 1.0, data(), nrows, o.data(), o.nrows, 0.0, out.data(), ncols);
//   return out;
// }

// //Matrix multiplication, right matrix is transposed
// matrixReal matrixReal::operator^(const matrixReal& o) const
// {
//   assert(ncols == o.ncols);
//   matrixReal out(nrows, o.nrows);
//   dgemm_("N","T", nrows, o.nrows, o.nrows, 1.0, data(), nrows, o.data(), o.nrows, 0.0, out.data(), nrows);
//   return out;
// }

matrixComp matrixComp:: operator+(const matrixComp& o) const
{
  assert(ncols == o.ncols && nrows == o.nrows);
  matrixComp out(nrows, o.ncols);
  transform(data(), data()+size(), o.data(), out.data(), plus<cplx>());
  return out;
}

matrixComp& matrixComp::operator+=(const matrixComp& o)
{
  *this = *this + o;
  return *this;
}

matrixComp matrixComp::operator-(const matrixComp& o) const
{
  assert(ncols == o.ncols && nrows == o.nrows);
  matrixComp out(nrows, o.ncols);
  transform(data(), data()+size(), o.data(), out.data(), minus<cplx>());
  return out;
}

matrixComp& matrixComp::operator-=(const matrixComp& o)
{
  *this = *this - o;
  return *this;
}

matrixComp matrixComp::operator*(const cplx& a) const
{
  matrixComp out(*this);
  out *= a;
  return out;
}

matrixComp& matrixComp::operator*=(const cplx& a)
{
  scale(a);
  return *this;
}

matrixComp matrixComp::operator/(const cplx& a) const
{
  matrixComp out(*this);
  out /= a;
  return out;
}

matrixComp& matrixComp::operator/=(const cplx& a)
{
  scale(1.0/a);
  return *this;
}

// double matrixReal::dot_product(const matrixReal& o) const {
//   assert(size() == o.size());
//   return ddot_(o.size(), data(), 1, o.data(), 1);
// }

// double matrixReal::norm() const {
//   return std::sqrt(dot_product(*this));
// }

// double matrixReal::rms() const {
//   return std::sqrt(dot_product(*this)/static_cast<double>(size()));
// }

// double matrixReal::variance() const {
//   return dot_product(*this)/static_cast<double>(size());
// }

void matrixComp::getEigvals(double* eigVals)
{
  assert (nrows == ncols);
  int info;
  int lwork = -1;
  std::unique_ptr <double[]> rwork (new double [nrows*3+2]);
  vector<cplx> wkopt (1,cplx(0.0));
  zheev_("N", "U", nrows, data(), nrows, eigVals, wkopt.data(), lwork, rwork.get(), info);
  lwork = int(real(wkopt[0]));
  std::unique_ptr <cplx[]> work (new cplx [lwork]);
  zheev_("N", "U", nrows, data(), nrows, eigVals, work.get(), lwork, rwork.get(), info);
  if (info > 0)
    throw std::runtime_error("Unable to diagonalize matrix");
  return;
}


// void matrixReal::diagonalize(double* eigVals, bool getLowEigVal, int keepNum)
// {
//   assert (nrows == ncols);
//   int info,iwork;
//   double work;

//   int lwork = -1;
//   int liwork = -1;
//   double abstol = 1.0e-8;
//   int il, iu;
//   matrixReal workEigVecs(nrows,keepNum);
//   vector<int> isuppz(std::max(1,keepNum),0);
//   if (getLowEigVal)
//   {
//       il = 1;
//       iu = keepNum;
//   }
//   else
//   {
//       il = nrows - keepNum + 1;
//       iu = nrows;
//   }
//   dsyevr_("V", "I", "U", nrows, data(), nrows, -20394857.0, 2345.2, il, iu, abstol, keepNum, eigVals, workEigVecs.data(), nrows, isuppz.data(), &work, lwork, &iwork, liwork, info);
//   lwork = work;

//   std::unique_ptr <double[]> workArray(new double [int(work)]);
//   liwork = iwork;
//   std::unique_ptr <int[]> iworkArray(new int [iwork]);

//   dsyevr_("V", "I", "U", nrows, data(), nrows, -20394857.0, 2345.2, il, iu, abstol, keepNum, eigVals, workEigVecs.data(), nrows, isuppz.data(), workArray.get(), lwork, iworkArray.get(), liwork, info);
//   if (info > 0)
//     throw std::runtime_error("Unable to diagonalize matrix with dyevr_");
//   zero();
//   setSub(0,0,workEigVecs);
//   return;
// }

// std::shared_ptr<matrixReal> matrixReal::transpose() const
// {
//   auto out = make_shared<matrixReal> (ncols,nrows);
//   mkl_domatcopy_("C","T",nrows,ncols,1.0,data(),nrows,out->data(),ncols);
//   return out;
// }

// tuple<shared_ptr<matrixReal>, shared_ptr<matrixReal>> matrixReal::svd(vector<double>& s)
// {
//   assert(s.size() >= std::min(ncols,nrows));
//   auto u = make_shared<matrixReal>(nrows, nrows);
//   auto vT = make_shared<matrixReal>(ncols, ncols);

//   int lwork = -1;
//   int info = 0;
//   dgesvd_("A","A", nrows, ncols, data(), nrows, s.data(), u->data(), nrows, vT->data(), ncols, s.data(), lwork, info);
//   lwork = s[0];
//   if (lwork <= 0)
//       throw runtime_error("dgesvd failed allocating lwork value");
//   unique_ptr<double[]> work(new double[lwork]);
//   dgesvd_("A","A", nrows, ncols, data(), nrows, s.data(), u->data(), nrows, vT->data(), ncols, work.get(), lwork, info);
//   if (info != 0)
//       throw runtime_error("dgesvd faied matrix decomposition");
//   return make_tuple(u,vT);
// }

// matrixReal matrixReal::kron(matrixReal &o) const
// {
//   matrixReal out(nrows*o.nrows,ncols*o.ncols);
//   for (int ii = 0; ii < nrows; ii++)
//   {
//     for (int jj = 0; jj < ncols; jj++)
//     {
//       out.setSub(ii*o.ncols,jj*o.nrows,o*element(ii,jj));
//     }
//   }
//   return out;
// }

// void matrixReal::ax_plus_y(const double a, matrixReal &o)
// {
//   assert(nrows*ncols == o.nrows*o.ncols);
//   daxpy_(size(), a, data(), 1, o.data(), 1);
// }
