#ifndef DMRG_UTILITIES
#define DMRG_UTILITIES

#include <complex>
#include <memory>
#include "mkl.h"

//BLAS
extern "C"
{
  void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
    const double* alpha, const double* a, const int* lda, const double* b, const int* ldb,
    const double* beta, double* c, const int* ldc);

  void dsyev_(const char*, const char*, const int*, double*, const int*, double*, double*, const int*, int*);

  double ddot_(const int*, const double*, const int*, const double*, const int*);

  void daxpy_(const int*, const double*, const double*, const int*, double*, const int*);
}

//LAPACK
extern "C"
{
 void dgesvd_(const char*, const char*, const int*, const int*, double*, const int*, double*, double*, const int*, double*, const int*,  double*, const int*, int*);

 //void dsyevr_(const char*, const char*, const char*, const int*, double*, const int*, const double*, const double*, const int*, const int*, const double*,  int*, double*, double*, const int*, int*, double*, int*, double*, int*, int*);
}

//DMRG interface, shamelessly lifted from BAGEL's src/util/f77.h file
namespace
{
  void dgemm_(const char* transa, const char* transb, const int m, const int n, const int k,
              const double alpha, const double* a, const int lda, const double* b, const int ldb,
              const double beta, double* c, const int ldc)
              { ::dgemm_(transa,transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc); }
  void dgemm_(const char* transa, const char* transb, const int m, const int n, const int k,
              const double alpha, const std::unique_ptr<double []>& a, const int lda, const std::unique_ptr<double []>& b, const int ldb,
              const double beta, std::unique_ptr<double []>& c, const int ldc)
             { ::dgemm_(transa,transb,&m,&n,&k,&alpha,a.get(),&lda,b.get(),&ldb,&beta,c.get(),&ldc); }

  void dsyev_(const char* a, const char* b, const int c, double* d, const int e, double* f, double* g, const int h, int& i)
             { ::dsyev_(a,b,&c,d,&e,f,g,&h,&i);}
  void dsyev_(const char* a, const char* b, const int c, std::unique_ptr<double []>& d, const int e,
             std::unique_ptr<double []>& f, std::unique_ptr<double []>& g, const int h, int& i)
             { ::dsyev_(a,b,&c,d.get(),&e,f.get(),g.get(),&h,&i);}

  void dgesvd_(const char* a, const char* b, const int c, const int d, double* e, const int f, double* g, double* h, const int i, double* j, const int k,
              double* l, const int m, int& n) { ::dgesvd_(a,b,&c,&d,e,&f,g,h,&i,j,&k,l,&m,&n); }


  void mkl_domatcopy_(const char* ordering, const char* trans, const int r, const int c, const double alpha,
                    const double* A, const int nr, double* B, const int nc)
                    {mkl_domatcopy(*ordering,*trans,r,c,alpha,A,nr,B,nc);}

  double ddot_(const int a, const double* b, const int c, const double* d, const int e) { return ::ddot_(&a,b,&c,d,&e);}

  double ddot_(const int a, const std::unique_ptr<double []>& b, const int c, const std::unique_ptr<double []>& d, const int e)
             { return ::ddot_(&a,b.get(),&c,d.get(),&e); }

  void dsyevr_(const char* jobz, const char* range, const char* uplo, const int n, double* a, const int lda, double vl,
                 double vu, const int il, const int iu, const double abstol, int m, double* w, double* z, const int ldz,
                int* isuppz, double* work, int lwork, int* iwork, int liwork, int& info)
                { ::dsyevr_(jobz, range, uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);}

  void daxpy_(const int a, const double b, const double* c, const int d, double* e, const int f) { ::daxpy_(&a,&b,c,&d,e,&f); }

}

#endif