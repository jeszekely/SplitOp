#ifndef SPLITOP_ARRAYS
#define SPLITOP_ARRAYS
#include <memory>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <random>
#include <complex>

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

//  Fill with a number
  void fill(T a) {std::fill_n(vals.get(), nx, a);}

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

//  Add a value to all elements
  template <typename U> void addVal(const U a)
  {
    std::for_each(data(), data()+size(), [&a](T& p){p+=a;});
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

// //Array-Array operations
  Array1D<T>& operator=(const Array1D<T>& o)
  {
    assert(nx == o.nx);
    std::copy_n(o.data(), o.size(), data());
    return *this;
  }

  template <typename U> Array1D operator*(const Array1D<U>& o) const
  {
  	assert(nx == o.nx);
  	Array1D<T> out(*this);
  	std::transform(out.data(), out.data()+out.size(), o.data(), out.data(), std::multiplies<T>());
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
  	std::transform(out.data(), out.data()+out.size(), o.data(), out.data(), std::divides<T>());
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
  	std::transform(out.data(), out.data()+out.size(), o.data(), out.data(), std::plus<T>());
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
  	std::transform(out.data(), out.data()+out.size(), o.data(), out.data(), std::minus<T>());
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
  		deriv = (1.0/(12.0*xstep))*(-25.0*vals[xx] + 48.0*vals[xx+1]-36.0*vals[xx+2]+16.0*vals[xx+3]-3.0*vals[xx+4]);
  	else if (xx > nx-4)
  		deriv = (1.0/(12.0*xstep))*(25.0*vals[xx]-48.0*vals[xx-1]+36.0*vals[xx-2]-16.0*vals[xx-3]+3.0*vals[xx-4]);
		else
			deriv = (1.0/(60.0*xstep))*(-1.0*vals[xx-3] + 9.0*vals[xx-2] - 45.0*vals[xx-1] + 45.0*vals[xx+1] - 9.0*vals[xx+2] + 1.0*vals[xx+3]);
		return deriv;
  }

  // T deriv_18(const int xx)

  // T deriv_24(const int xx)

  //2nd derivative, 4th order accuracy
  T deriv_24(const int xx)
  {
  	T deriv;
  	if (xx < 3)
			deriv = (1.0/(12.0*xstep*xstep))*(45.0*vals[xx] - 154.0*vals[xx+1] + 214.0*vals[xx+2] - 156.0*vals[xx+3] + 61.0*vals[xx+4] - 10.0*vals[xx+5]);
		else if (xx > nx-4)
			deriv = (1.0/(12.0*xstep*xstep))*(45.0*vals[xx] - 154.0*vals[xx-1] + 214.0*vals[xx-2] - 156.0*vals[xx-3] + 61.0*vals[xx-4] - 10.0*vals[xx-5]);
		else
			deriv = (1.0/(180.0*xstep*xstep))*(2.0*vals[xx-3] - 27.0*vals[xx-2] + 270.0*vals[xx-1] - 490.0*vals[xx] + 270.0*vals[xx+1] - 27.0*vals[xx+2] + 2.0*vals[xx+3]);
		return deriv;
	}
  // T deriv_28(const int xx)

  T integrate_rect()
  {
  	return xstep*std::accumulate(data(),data()+size(),T(0.0));
  }

  T integrate_rect(const int xL, const int xU)
  {
    return xstep*std::accumulate(data()+xL,data()+xL+xU,T(0.0));
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
			out << std::scientific << std::setprecision(3) << o(elem) << "\t";
	}
	else
	{
		for (int elem = 0; elem < 5; elem++)
			out <<  std::scientific << std::setprecision(3) << o(elem) << "\t";
		out << "...\t";
		for (int elem = o.nx-5; elem < o.nx; elem++)
			out <<  std::scientific << std::setprecision(3) << o(elem) << "\t";
	}
	out << std::endl;
  return out;
};

template <class T>
void print1D(int nn, std::ostream &out,T& o)
{
	out << o(nn) << "\t";
};

template <class A, class... B>
void print1D(int nn, std::ostream &out, A head, B... tail)
{
	print1D(nn,out,head);
	print1D(nn,out,tail...);
	out << std::endl;
};

template <class A, class... B>
void printArrays(int nn, std::ostream &out, A head, B... tail)
{
	for (int ii = 0; ii<nn; ii++)
	{
		print1D(ii,out,head);
		print1D(ii,out,tail...);
		out << std::endl;
	}
};

template <typename T> unsigned int Array1D<T>::memSize = 0;


/***********************************

Begin 2D Array class definitions

************************************/

template <typename T> class Array2D
{
protected:
  static unsigned int memSize;
  size_t nx;
  size_t ny;
  T xstep, ystep;
  std::unique_ptr<T[]> vals;

public:
//  Constructor
  Array2D(const int NX, const int NY, T xs, T ys) : nx(NX), ny(NY), xstep(xs), ystep(ys), vals(std::unique_ptr<T[]>(new T[NX*NY]))
  {
    zero();
    memSize += sizeof(T)*size(); //number of bytes allocated to instance of class
  }
//  Copy Constructor
  Array2D(const Array2D& o) : nx(o.nx), ny(o.ny), xstep(o.xstep), ystep(o.ystep), vals(std::unique_ptr<T[]>(new T[nx*ny]))
  {
    std::copy_n(o.vals.get(), nx*ny, vals.get());
    memSize += sizeof(T)*size(); //number of bytes allocated to instance of class
  }
  //  Move Constructor
  Array2D(Array2D&& o) : nx(o.nx), ny(o.ny), xstep(o.xstep), ystep(o.ystep), vals(std::move(o.vals)) {o.nx = 0; o.ny = 0;};

//  Destructor
  ~Array2D(){memSize -= sizeof(T)*size();} //number of bytes allocated to instance of class

//  Access functions
  size_t size() const { return nx*ny; }
  T* data() { return vals.get(); }
  const T* data() const { return vals.get(); }

//  Fill with zeroes
  void zero() {std::fill_n(vals.get(), nx*ny, T(0.0));}

//  Fill with a number
  void fill(T a) {std::fill_n(vals.get(), nx*ny, a);}

  size_t Nx() const {return nx;}
  size_t Ny() const {return ny;}

  void random()
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> dis(-2, 2);
    std::generate_n(vals.get(), nx*ny, [&dis, &gen](){return dis(gen);});
  }

//  Accessor functions
  T& element(const int xx, const int yy)
  {
    return vals[xx + yy*nx];
  }

  const T& element(const int xx, const int yy) const
  {
    return vals[xx + yy*nx];
  }

  T& operator()(const int xx, const int yy)
  {
    return element(xx,yy);
  }

  const T& operator()(const int xx, const int yy) const
  {
    return element(xx,yy);
  }

//  Apply a scaling factor to all elements
  template <typename U> void scale(const U a)
  {
    std::for_each(data(), data()+size(), [&a](T& p){p*=a;});
  }

//  Add a value to all elements
  template <typename U> void addVal(const U a)
  {
    std::for_each(data(), data()+size(), [&a](T& p){p+=a;});
  }

//  Print memory usage for all matrices
  void printMem() const
  {
    std::cout << "Current memory allocated to this array: " << size()*sizeof(T) << " bytes." << std::endl;
    std::cout << "Total memory allocated for data storage: " << memSize << " bytes." << std::endl;
    return;
  }

//Extract a portion of the matrix starting at element(r,c), get (nr x nc matrix)
  Array2D <T> getSub(const int xindex, const int yindex, const int nxvals, const int nyvals) const
  {
    assert(xindex + nxvals <= nx && yindex + nyvals <= ny);
    Array2D <T> out(nxvals,nyvals,xstep,ystep);
    for (int ii = 0; ii < nyvals; ii++)
      std::copy_n(&element(xindex,yindex+ii),nxvals,&out.element(0,ii));
    return out;
 }

//Place matrix o at position (r,c)
  template <class U>
  void setSub(int xindex, int yindex, U o)
  {
    assert(xindex + o.nx <= nx && yindex + o.ny <= ny);
    for (int ii = 0; ii < o.ny; ii++)
      std::copy_n(&o(0,ii), o.nx, &element(xindex,yindex+ii));
    return;
  }

  Array2D operator=(const Array2D<T>& o)
  {
    assert(nx == o.nx && ny = o.ny);
    std::copy_n(o.data(), o.size(), data());
    return *this;
  }

//element by element operations
  template <typename U> Array2D operator*(const Array2D<U>& o) const
  {
    assert(nx == o.nx && ny == o.ny);
    Array2D<T> out(*this);
    std::transform(out.data(), out.data()+out.size(), o.data(), out.data(), std::multiplies<T>());
    return out;
  }

  template <typename U> Array2D& operator*=(const Array2D<U>& o)
  {
    *this = *this * o;
    return *this;
  }

  template <typename U> Array2D operator/(const Array2D<U>& o) const
  {
    assert(nx == o.nx && ny == o.ny);
    Array2D<T> out(*this);
    std::transform(out.data(), out.data()+out.size(), o.data(), out.data(), std::divides<T>());
    return out;
  }
  template <typename U> Array2D& operator/=(const Array2D<U>& o)
  {
    *this = *this / o;
    return *this;
  }

  template <typename U> Array2D operator+(const Array2D<U>& o) const
  {
    assert(nx == o.nx && ny == o.ny);
    Array2D<T> out(*this);
    std::transform(out.data(), out.data()+out.size(), o.data(), out.data(), std::plus<T>());
    return out;
  }
  template <typename U> Array2D& operator+=(const Array2D<U>& o)
  {
    *this = *this + o;
    return *this;
  }

  template <typename U> Array2D operator-(const Array2D<U>& o) const
  {
    assert(nx == o.nx && ny == o.ny);
    Array2D<T> out(*this);
    std::transform(out.data(), out.data()+out.size(), o.data(), out.data(), std::minus<T>());
    return out;
  }

  template <typename U> Array2D& operator-=(const Array2D<U>& o)
  {
    *this = *this - o;
    return *this;
  }

//Scalar-Array Operations
  template <typename U> Array2D operator*(const U& a) const
  {
    Array2D<T> out(*this);
    out *= a;
    return out;
  }

  template <typename U> Array2D operator/(const U& a) const
  {
    Array2D<T> out(*this);
    out /= a;
    return out;
  }
  template <typename U> void operator*=(const U& a) {scale(a);}
  template <typename U> void operator/=(const U& a) {scale(1/a);}

  // T deriv_14(const int xx)

  // 1st derivative, 6th order accuracy, x direction
  T deriv_16_x(const int xx, const int yy)
  {
    T deriv;
    if (xx < 3)
      deriv = (1.0/(12.0*xstep))*(-25.0*element(xx,yy) + 48.0*element(xx+1,yy)-36.0*element(xx+2,yy)+16.0*element(xx+3,yy)-3.0*element(xx+4,yy));
    else if (xx > nx-4)
      deriv = (1.0/(12.0*xstep))*(25.0*element(xx,yy)-48.0*element(xx-1,yy)+36.0*element(xx-2,yy)-16.0*element(xx-3,yy)+3.0*element(xx-4,yy));
    else
      deriv = (1.0/(60.0*xstep))*(-1.0*element(xx-3,yy) + 9.0*element(xx-2,yy) - 45.0*element(xx-1,yy) + 45.0*element(xx+1,yy) - 9.0*element(xx+2,yy) + 1.0*element(xx+3,yy));
    return deriv;
  }

  // 1st derivative, 6th order accuracy, y direction
  T deriv_16_y(const int xx, const int yy)
  {
    T deriv;
    if (yy < 3)
      deriv = (1.0/(12.0*ystep))*(-25.0*element(xx,yy) + 48.0*element(xx,yy+1)-36.0*element(xx,yy+2)+16.0*element(xx,yy+3)-3.0*element(xx,yy+4));
    else if (yy > ny-4)
      deriv = (1.0/(12.0*ystep))*(25.0*element(xx,yy)-48.0*element(xx,yy-1)+36.0*element(xx,yy-2)-16.0*element(xx,yy-3)+3.0*element(xx,yy-4));
    else
      deriv = (1.0/(60.0*ystep))*(-1.0*element(xx,yy-3) + 9.0*element(xx,yy-2) - 45.0*element(xx,yy-1) + 45.0*element(xx,yy+1) - 9.0*element(xx,yy+2) + 1.0*element(xx,yy+3));
    return deriv;
  }

  // T deriv_18(const int xx)

  // T deriv_24(const int xx)

  //2nd derivative, 4th order accuracy
  T deriv_24_xx(const int xx, const int yy)
  {
    T deriv;
    if (xx < 3)
      deriv = (1.0/(12.0*xstep*xstep))*(45.0*element(xx,yy) - 154.0*element(xx+1,yy) + 214.0*element(xx+2,yy) - 156.0*element(xx+3,yy) + 61.0*element(xx+4,yy) - 10.0*element(xx+5,yy));
    else if (xx > nx-4)
      deriv = (1.0/(12.0*xstep*xstep))*(45.0*element(xx,yy) - 154.0*element(xx-1,yy) + 214.0*element(xx-2,yy) - 156.0*element(xx-3,yy) + 61.0*element(xx-4,yy) - 10.0*element(xx-5,yy));
    else
      deriv = (1.0/(180.0*xstep*xstep))*(2.0*element(xx-3,yy) - 27.0*element(xx-2,yy) + 270.0*element(xx-1,yy) - 490.0*element(xx,yy) + 270.0*element(xx+1,yy) - 27.0*element(xx+2,yy) + 2.0*element(xx+3,yy));
    return deriv;
  }

  T deriv_24_yy(const int xx, const int yy)
  {
    T deriv;
    if (yy < 3)
      deriv = (1.0/(12.0*ystep*ystep))*(45.0*element(xx,yy) - 154.0*element(xx,yy+1) + 214.0*element(xx,yy+2) - 156.0*element(xx,yy+3) + 61.0*element(xx,yy+4) - 10.0*element(xx,yy+5));
    else if (yy > ny-4)
      deriv = (1.0/(12.0*ystep*ystep))*(45.0*element(xx,yy) - 154.0*element(xx,yy-1) + 214.0*element(xx,yy-2) - 156.0*element(xx,yy-3) + 61.0*element(xx,yy-4) - 10.0*element(xx,yy-5));
    else
      deriv = (1.0/(180.0*ystep*ystep))*(2.0*element(xx,yy-3) - 27.0*element(xx,yy-2) + 270.0*element(xx,yy-1) - 490.0*element(xx,yy) + 270.0*element(xx,yy+1) - 27.0*element(xx,yy+2) + 2.0*element(xx,yy+3));
    return deriv;
  }
  // T deriv_28(const int xx)


  T integrate_rect()
  {
    return xstep*ystep*std::accumulate(data(),data()+size(),T(0.0));
  }

  T integrate_rect(const int xL, const int xU, const int yL, const int yU)
  {
    T sum = 0;
    for (int ii = yL; ii < yU; ii++)
      sum += std::accumulate(&element(xL,ii),&element(xL+xU,ii),T(0.0));
    return xstep*ystep*sum;
  }

  // T integrate_trap()
  // T integrate_simp()

  template <typename U> friend std::ostream &operator<<(std::ostream &out, const Array2D <U> &o);

};

//Overload the << operator to print a matrix
template <typename T>
std::ostream &operator<<(std::ostream &out, const Array2D <T> &o)
{
  for (int row = 0; row < std::min(10,int(o.nx)); row++)
  {
    for (int col = 0; col < std::min(10,int(o.ny)); col++)
    {
      out << std::setprecision(3) << o(row,col) << "\t";
    }
    out << "\n";
  }
  out << std::endl;
  return out;
};

template <typename T> unsigned int Array2D<T>::memSize = 0;
#endif
