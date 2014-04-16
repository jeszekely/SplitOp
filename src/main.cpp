#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <ctime>

#include "wvfxn.hpp"
#include "input_parser.hpp"

using namespace std;
int main(int argc, char const *argv[])
{
  programInputs IP("inputs.json");
  cout << IP.nx << endl;

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
  cout << endl;
  Array2D <float> AA(20,20,0.1,0.1);
  AA.zero();
  cout << AA;
  AA.fill(1.0);
  cout << AA << AA.getSub(0,0,3,3);
  cout << AA.integrate_rect() << endl;

  return(0);
}
