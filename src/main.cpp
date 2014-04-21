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


#if 0
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
#endif

//Set up basics for the 2D calculation
SplitOp2D MainCalc(IP);
MainCalc.wvfxn->mass1 = IP.m_electron;
MainCalc.wvfxn->mass2 = IP.m_C60;
MainCalc.initializeTDSE([&](double a, double b){return HPotential(a,b,IP);}, [&](double a, double b){return HKinetic2D(a,b,IP.m_electron,IP.m_C60);});

//Setting up electronic subspace calculation
int JL = -1;
int JR = -1;
for (int ii = 0; ii < IP.nx; ii++)
{
  if (JL < 0 && MainCalc.xgrid->element(ii) >= IP.requil - 45.0)
    JL = ii;
  if (JR < 0 && MainCalc.xgrid->element(ii) >= IP.requil + 45.0)
    JR = ii;
  if (JR > 0 && JL > 0)
    break;
}
int nxJunction = JR - JL;
if (nxJunction % 2 == 1)
{
  nxJunction++;
  JL++;
}
programInputs JuncIP = IP;
JuncIP.nx = nxJunction;
JuncIP.xmin = MainCalc.xgrid->element(JL);
JuncIP.xmax = MainCalc.xgrid->element(JR);
cout << "The junction spans from index " << JL << " to " << JR << " in the electronic subspace." << endl;
SplitOp1D ElecCalc(JuncIP);
ElecCalc.initializeITP([&](double a){return ElecParab(a,IP.requil,IP);}, [&](double a){return HKinetic1D(a,IP.m_electron);});
Array2D <cplx> ElecStates(nxJunction,IP.nelec,JuncIP.xmin,ElecCalc.xstep);
ElecCalc.propagateITP(ElecStates);
if (IP.elecstates == 1)
{
  ofstream ElecITP;
  ElecITP.open("output_data/ElecITP_Results.txt");
  ElecITP << "Junction Coordinate (au) \t Potential \t States" << endl;
  for (int kk = 0; kk < ElecStates.Nx(); kk++)
  {
    ElecITP << ElecCalc.xgrid->element(kk) << "\t" << real(ElecCalc.Vgrid->element(kk));
    for (int jj = 0; jj < ElecStates.Ny(); jj++)
      ElecITP << "\t" << real(ElecStates(kk,jj));
    ElecITP << endl;
  }
  ElecITP.close();
}


  return(0);
}
