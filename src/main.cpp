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
#if 0
  SplitOp1D TestCalc(IP);
  transform(TestCalc.xgrid->data(),TestCalc.xgrid->data()+IP.nx,TestCalc.wvfxn->data(),[&](double x){return cplx(wvfxnElectron(x,IP));});
  TestCalc.wvfxn->normalize();

  TestCalc.initializeTDSE([&](double a){return cplx(0.0,-1.0)*absorbingPotential(a,IP)+Ve(a,IP);}, [&](double a){return HKinetic1D(a,IP.m_electron);});
  TestCalc.propagateStep(100);
  printArrays(7000,cout,*TestCalc.xgrid,*TestCalc.wvfxn);
#endif
  //ITP test
  IP.xmin = -10.0;
  IP.xmax = 10.0;
  IP.xi = 2;
  IP.elecpos = -2.0;
  SplitOp1D ITPCalc(IP);
  transform(ITPCalc.xgrid->data(),ITPCalc.xgrid->data()+IP.nx,ITPCalc.wvfxn->data(),[&](double x){return cplx(wvfxnElectron(x,IP));});
  ITPCalc.wvfxn->normalize();
  ITPCalc.initializeITP([&](double a){return IP.m_electron*pow(a,2);}, [&](double a){return HKinetic1D(a,IP.m_electron);});
  Array2D <cplx> States(IP.nx,3,0.0,0.0);
  ITPCalc.propagateITP(States,1.0e-4);

  wvfxn1D GS(*ITPCalc.wvfxn);
  wvfxn1D E1(*ITPCalc.wvfxn);
  wvfxn1D E2(*ITPCalc.wvfxn);

  copy_n(&States(0,0),IP.nx,GS.data());
  copy_n(&States(0,1),IP.nx,E1.data());
  copy_n(&States(0,2),IP.nx,E2.data());

  GS.normalize();
  E1.normalize();
  E2.normalize();

  printArrays(7000,cout,*ITPCalc.xgrid,GS,E1,E2);

  return(0);
}
