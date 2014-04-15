#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <ctime>

#include "arrays.hpp"

using namespace std;
int main(int argc, char const *argv[])
{
  Array1D <double> A(100,0,0);
  for (int ii = 0; ii<A.Nx(); ii++)
    A(ii) = -10.0;
  cout << A.integrate_rect() << endl;
  Array1D  <double> X(100,0,0.5);
  X.fill_array();
  cout << A << X;
  Array1D <double> Y(X);
  for (int ii = 0; ii<X.Nx(); ii++)
    Y(ii) = pow(X(ii),2);
  cout << X(25) << "\t" << Y.deriv_16(25) << endl;
  Y*=10.0;
  cout << X(25) << "\t" << Y.deriv_16(25) << endl;
  A.random();
  cout << A;
  Array1D <double> C(10,0,0);
  C.random();
  cout << C << endl;
  A.setSub(2,C);
  cout << A << endl << A.getSub(0,4);

  printArrays(10,cout,A,X);
  return(0);
}
