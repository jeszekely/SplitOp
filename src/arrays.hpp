#ifndef SPLITOP_ARRAYS
#define SPLITOP_ARRAYS

#include <memory>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <random>

template <typename T> class Array1D
{
protected:
  static unsigned int memSize;
  size_t nx;
  T xinit, xstep;
  std::unique_ptr<T[]> vals;

public:
//  Constructor
  Array1D(const int NX, T xi, T xs) : nx(NX), xinit(xi), xstep(xs), vals(std::unique_ptr<T[]>(new T[NX]))
  {
    zero();
    memSize += sizeof(T)*size(); //number of bytes allocated to instance of class
  }
//  Copy Constructor
  Array1D(const Array1D& o) : nx(o.nx), xinit(o.xinit), xstep(o.xstep), vals(std::unique_ptr<T[]>(new T[nx]))
  {
    std::copy_n(o.vals.get(), nx, vals.get());
    memSize += sizeof(T)*size(); //number of bytes allocated to instance of class
  }
  //  Move Constructor
  Array1D(Array1D&& o) : nx(o.nx), xinit(o.xinit), xstep(o.xstep), vals(std::move(o.vals)) {o.nx = 0;};

//  Destructor
  ~Array1D(){memSize -= sizeof(T)*size();} //number of bytes allocated to instance of class

//  Access functions
  size_t size() const { return nx; }
  T* data() { return vals.get(); }
  const T* data() const { return vals.get(); }

//  Fill with zeroes
  void zero() {std::fill_n(vals.get(), nx, T(0.0));}

	size_t Nx() const {return nx;}

  void random()
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> dis(-2, 2);
    std::generate_n(vals.get(), nx, [&dis, &gen](){return dis(gen);});
  }

  void fill_array()
  {
  	for (int ii = 0; ii < nx; ii++)
  		element(ii) = ii*xstep + xinit;
  	return;
  }

//  Accessor functions
  T& element(const int nn)
  {
    return vals[nn];
  }

  const T& element(const int nn) const
  {
    return vals[nn];
  }

  T& operator()(const int nn)
  {
    return element(nn);
  }

  const T& operator()(const int nn) const
  {
    return vals[nn];
  }

//  Apply a scaling factor to all elements
  template <typename U> void scale(const U a)
  {
    std::for_each(data(), data()+size(), [&a](T& p){p*=a;});
  }

//  Print memory usage for all matrices
  void printMem() const
  {
    std::cout << "Current memory allocated to this array: " << size()*sizeof(T) << " bytes." << std::endl;
    std::cout << "Total memory allocated for data storage: " << memSize << " bytes." << std::endl;
    return;
  }

//Extract a portion of the matrix starting at element(r,c), get (nr x nc matrix)
  Array1D <T> getSub(const int index, const int nvals) const
  {
    assert(index + nvals <= nx);
    Array1D <T> out(nvals,xinit,xstep);
    std::copy_n(&element(index),nvals,&out.element(0));
    return out;
 }

//Place matrix o at position (r,c)
  template <class U>
  void setSub(int index, U o)
  {
    assert(index + o.nx <= nx);
    std::copy_n(&o(0), o.nx, &element(index));
    return;
  }

//Array-Array operations
  template <typename U> Array1D operator*(const Array1D<U>& o) const
  {
  	assert(nx == o.nx);
  	Array1D<T> out(*this);
  	std::transform(out.data(), out.data()+out.size(), o.data(), std::multiplies<T>());
  	return out;
  }
  template <typename U> Array1D& operator*=(const Array1D<U>& o)
  {
  	*this = *this * o;
  	return *this;
  }

  template <typename U> Array1D operator/(const Array1D<U>& o) const
  {
  	assert(nx == o.nx);
  	Array1D<T> out(*this);
  	std::transform(out.data(), out.data()+out.size(), o.data(), std::divides<T>());
  	return out;
  }
  template <typename U> Array1D& operator/=(const Array1D<U>& o)
  {
  	*this = *this / o;
  	return *this;
  }

  template <typename U> Array1D operator+(const Array1D<U>& o) const
  {
  	assert(nx == o.nx);
  	Array1D<T> out(*this);
  	std::transform(out.data(), out.data()+out.size(), o.data(), std::plus<T>());
  	return out;
  }
  template <typename U> Array1D& operator+=(const Array1D<U>& o)
  {
  	*this = *this + o;
  	return *this;
  }

  template <typename U> Array1D operator-(const Array1D<U>& o) const
  {
  	assert(nx == o.nx);
  	Array1D<T> out(*this);
  	std::transform(out.data(), out.data()+out.size(), o.data(), std::minus<T>());
  	return out;
  }
  template <typename U> Array1D& operator-=(const Array1D<U>& o)
  {
  	*this = *this - o;
  	return *this;
  }

//Scalar-Array Operations
  template <typename U> Array1D operator*(const U& a) const
  {
  	Array1D<T> out(*this);
  	out *= a;
  	return out;
  }

  template <typename U> Array1D operator/(const U& a) const
  {
  	Array1D<T> out(*this);
  	out /= a;
  	return out;
  }
  template <typename U> void operator*=(const U& a) {scale(a);}
  template <typename U> void operator/=(const U& a) {scale(1/a);}

  // T deriv_14(const int xx)

  //1st derivative, 6th order accuracy
  T deriv_16(const int xx)
  {
  	T deriv;
  	if (xx < 3)
  		deriv = (1/(12*xstep))*(-25*vals[xx] + 48*vals[xx+1]-36*vals[xx+2]+16*vals[xx+3]-3*vals[xx+4]);
  	else if (xx > nx-4)
  		deriv = (1/(12.0*xstep))*(25.0*vals[xx]-48.0*vals[xx-1]+36.0*vals[xx-2]-16.0*vals[xx-3]+3.0*vals[xx-4]);
		else
			deriv = (1.0/(60*xstep))*(-1.0*vals[xx-3] + 9.0*vals[xx-2] - 45.0*vals[xx-1] + 45.0*vals[xx+1] - 9.0*vals[xx+2] + 1.0*vals[xx+3]);
		return deriv;
  }

  // T deriv_18(const int xx)

  // T deriv_24(const int xx)

  //2nd derivative, 6th order accuracy
  T deriv_26(const int xx)
  {
  	T deriv;
  	if (xx < 3)
			deriv = (1.0/(12*xstep*xstep))*(45*vals[xx] - 154*vals[xx+1] + 214*vals[xx+2] - 156*vals[xx+3] + 61*vals[xx+4] - 10*vals[xx+5]);
		else if (xx > nx-4)
			deriv = (1.0/(12*xstep*xstep))*(45*vals[xx] - 154*vals[xx-1] + 214*vals[xx-2] - 156*vals[xx-3] + 61*vals[xx-4] - 10*vals[xx-5]);
		else
			deriv = (1.0/(180*xstep*xstep))*(2*vals[xx-3] - 27*vals[xx-2] + 270*vals[xx-1] - 490*vals[xx] + 270*vals[xx+1] - 27*vals[xx+2] + 2*vals[xx+3]);
		return deriv;
	}
  // T deriv_28(const int xx)

  T integrate_rect()
  {
  	return std::accumulate(data(),data()+size(),0.0);
  }
  // T integrate_trap()
  // T integrate_simp()

  template <typename U> friend std::ostream &operator<<(std::ostream &out, const Array1D <U> &o);

  template <typename... Args> friend void print1D(std::ostream &out, const Args&... o);
};

//Overload the << operator to print a matrix
template <typename T>
std::ostream &operator<<(std::ostream &out, const Array1D <T> &o)
{
	if (o.nx < 10)
	{
		for (int elem = 0; elem < o.nx; elem++)
			out << std::setprecision(3) << o(elem) << "\t";
	}
	else
	{
		for (int elem = 0; elem < 5; elem++)
			out << std::setprecision(3) << o(elem) << "\t";
		out << "...\t";
		for (int elem = o.nx-5; elem < o.nx; elem++)
			out << std::setprecision(3) << o(elem) << "\t";
	}
	out << std::endl;
  return out;
};

template <class T>
void print1D(int nn, std::ostream &out,T& o)
{
	out << o(nn) << "\t";
}

template <class A, class... B>
void print1D(int nn, std::ostream &out, A head, B... tail)
{
	print1D(nn,out,head);
	print1D(nn,out,tail...);
	out << std::endl;
}

template <class A, class... B>
void printArrays(int nn, std::ostream &out, A head, B... tail)
{
	for (int ii = 0; ii<nn; ii++)
	{
		print1D(ii,out,head);
		print1D(ii,out,tail...);
		out << std::endl;
	}
}



template <typename T> unsigned int Array1D<T>::memSize = 0;
#endif
