#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <ctime>

#include "wvfxn.hpp"
#include "input_parser.hpp"
#include "splitop.hpp"
#include "junction.hpp"

using namespace std;
int main(int argc, char const *argv[])
{
  programInputs IP("inputs.json");
  SplitOp1D TestCalc(IP);
  transform(TestCalc.xgrid->data(),TestCalc.xgrid->data()+IP.nx,TestCalc.wvfxn->data(),[&](double x){return cplx(wvfxnElectron(x,IP));});
  TestCalc.wvfxn->normalize();

  TestCalc.initializeTDSE([&](double a){return cplx(0.0,-1.0)*absorbingPotential(a,IP)+Ve(a,IP);}, [&](double a){return HKinetic1D(a,IP.m_electron);});
  TestCalc.propagateNSteps(100);
  printArrays(7000,cout,*TestCalc.xgrid,*TestCalc.wvfxn);

  return(0);
}
