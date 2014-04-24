#ifndef SPLITOP_POLYNOMIAL
#define SPLITOP_POLYNOMIAL

#include <assert.h>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>

template <typename T>
class polynomial
{
public:
  std::shared_ptr<std::vector<T>> vals;

  polynomial(const int n) : vals(std::make_shared<std::vector<T>>(n))
  {
    fill_n(vals->data(),size(),T(0.0));
  }
  polynomial(const polynomial& o) : vals(std::make_shared<std::vector<T>>(o.size()))
  {
    std::copy_n(o.vals->data(), size(), vals->data());
  }
  // polynomial(polynomial&& o) : size()(o.size()), vals(std::move(o.vals)){}

  T& element(const int nn)
  {
    return vals->at(nn);
  }

  const T& element(const int nn) const
  {
    return vals->at(nn);
  }

  T& operator()(const int nn)
  {
    return element(nn);
  }

  const T& operator()(const int nn) const
  {
    return element(nn);
  }

  const int size() const { return vals->size(); }

  polynomial& operator=(const polynomial<T>& o)
  {
    vals = std::make_shared<std::vector<T>>(o.size());
    copy_n(&o(0),o.size(),&element(0));
    return *this;
  }

  polynomial operator+(const polynomial<T>& o) const
  {
    int len;
    const polynomial<T> *p;
    size() >= o.size() ? (p=this, len=size()) : (p=&o, len=o.size());
    polynomial<T> out(len);
    std::transform(vals->data(),vals->data()+std::min(size(),o.size()),o.vals->data(),out.vals->data(),[](T a, T b){return a+b;});
    if (size() != o.size()) std::copy_n(&p->element(std::min(size(),o.size())),abs(size()-o.size()),&out(std::min(size(),o.size())));
    return out;
  }

  polynomial& operator+=(const polynomial<T>& o)
  {
    *this = *this + o;
    return *this;
  }

  polynomial operator-(const polynomial<T>& o) const
  {
    int len, sign;
    int min = std::min(size(),o.size());
    const polynomial<T> *p;
    size() >= o.size() ? (p=this, len=size(), sign=1) : (p=&o, len=o.size(), sign=-1);
    polynomial<T> out(len);
    std::transform(vals->data(),vals->data()+min,o.vals->data(),out.vals->data(),[](T a, T b){return a-b;});
    if (size() != o.size())
    {
      std::copy_n(&p->element(min),abs(size()-o.size()),&out(min));
      if (sign == -1) std::transform(&out(min),&out(min)+abs(size()-o.size()),&out(min),[](T a){return -1.0*a;});
    }
    return out;
  }

  polynomial& operator-=(const polynomial<T>& o)
  {
    *this = *this - o;
    return *this;
  }

  polynomial operator*(const polynomial<T>& o) const
  {
    int len = size()+o.size()-1;
    polynomial<T> out(len);
    for (int ii = 0; ii<size(); ii++)
    {
      for (int jj = 0; jj<o.size(); jj++)
        out(ii+jj) += (element(ii)*o(jj));
    }
    return out;
  }

  polynomial& operator*=(const polynomial<T>& o)
  {
    *this = *this * o;
    return *this;
  }

  template <typename U>
  void scale(const U a)
  {
    std::for_each(vals->data(), vals->data()+size(), [&a](T& p){p*=a;});
    return;
  }

};

template <typename T>
std::ostream &operator<<(std::ostream &out, const polynomial<T> &o)
{
  for (int ii = 0; ii < o.size(); ii++)
    out << o(ii) << "\t";
  out << endl;
  return out;
};

//"probabilist's" Hermite polynomials
template <typename T>
std::shared_ptr<polynomial<T>> Hermite(int nn)
{
  std::shared_ptr<polynomial<T>> p0,p1,pn;
  polynomial<T> a1(2);
  a1(1)          = 1;
  p0             = make_shared<polynomial<T>>(1);
  p0->element(0) = T(1.0);
  p1             = make_shared<polynomial<T>>(2);
  p1->element(1) = T(1.0);
  if (nn == 0) return p0;
  if (nn == 1) return p1;
  int kk = 1;
  while (kk < nn)
  {
    kk++;
    p0->scale(kk-1);
    pn = make_shared<polynomial<T>>(*p1*a1-(*p0));
    p0 = p1;
    p1 = pn;
  }
  return pn;
};

//Chebyshev polynomials of the first kind
template <typename T>
std::shared_ptr<polynomial<T>> Chebyshev(int nn)
{
  std::shared_ptr<polynomial<T>> p0,p1,pn;
  polynomial<T> a1(2);
  a1(1)          = 2;
  p0             = make_shared<polynomial<T>>(1);
  p0->element(0) = T(1.0);
  p1             = make_shared<polynomial<T>>(2);
  p1->element(1) = T(1.0);
  if (nn == 0) return p0;
  if (nn == 1) return p1;
  int kk = 1;
  while (kk < nn)
  {
    kk++;
    pn = make_shared<polynomial<T>>(*p1*a1-(*p0));
    p0 = p1;
    p1 = pn;
  }
  return pn;
};



#endif