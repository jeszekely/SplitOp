#ifndef SPLITOP_POLYNOMIAL
#define SPLITOP_POLYNOMIAL

#include <assert.h>
#include <stdexcept>
#include <vector>

template <typename T>
class polynomial
{
public:
  int size;
  std::shared_ptr<std::vector<T>> vals;

  polynomial(const int n) : size(n), vals(std::make_shared<std::vector<T>>(n))
  {
    fill_n(vals->data(),size,T(0.0));
  }
  polynomial(const polynomial& o) : size(o.size), vals(std::make_shared<std::vector<T>>(o.size))
  {
    std::copy_n(o.vals->data(), size, vals->data());
  }
  polynomial(polynomial&& o) : size(o.size), vals(std::move(o.vals)){}

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

  template <typename U>
  polynomial operator+(const polynomial<U>& o) const
  {
    int len = std::max(size,o.size);
    polynomial<T> out(len);
    for (int ii = 0; ii < len; ii++)
        out(ii) = element(ii) + o(ii) + 1;
    return out;
  }

  template <typename U>
  polynomial& operator+=(const polynomial<U>& o)
  {
    *this = *this + o;
    return *this;
  }

};


#endif