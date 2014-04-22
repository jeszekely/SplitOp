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

//Set up basics for the 2D calculation
  SplitOp2D MainCalc(IP);
  MainCalc.wvfxn->mass1 = IP.m_electron;
  MainCalc.wvfxn->mass2 = IP.m_C60;
  MainCalc.initializeTDSE([&](double a, double b){return HPotential(a,b,IP);}, [&](double a, double b){return HKinetic2D(a,b,IP.m_electron,IP.m_C60);});

/*******************************************
Electronic Supspace Calculations
********************************************/

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
  JuncIP.nx            = nxJunction;
  JuncIP.xmin          = MainCalc.xgrid->element(JL);
  JuncIP.xmax          = MainCalc.xgrid->element(JR);
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

/*******************************************
Molecular/Phonon Subspace Calculations
********************************************/
  programInputs MoleIP = JuncIP;
  MoleIP.nx            = IP.ny;
  MoleIP.xmin          = IP.ymin;
  MoleIP.ymax          = IP.ymax;

  SplitOp1D MoleCalc(MoleIP);
  MoleCalc.initializeITP([&](double a){return VM(a,IP);}, [&](double a){return HKinetic1D(a,IP.m_C60);});
  Array2D <cplx> MoleStates(IP.ny,IP.nphonon,IP.ymin,MoleCalc.xstep);
  MoleCalc.propagateITP(MoleStates);
  if (IP.phononstates == 1)
  {
    ofstream PhononITP;
    PhononITP.open("output_data/PhononITP_Results.txt");
    PhononITP << "Phonon Coordinate (au) \t Potential \t States" << endl;
    for (int kk = 0; kk < MoleStates.Nx(); kk++)
    {
      PhononITP << MoleCalc.xgrid->element(kk) << "\t" << real(MoleCalc.Vgrid->element(kk));
      for (int jj = 0; jj < MoleStates.Ny(); jj++)
        PhononITP << "\t" << real(MoleStates(kk,jj));
      PhononITP << endl;
    }
    PhononITP.close();
  }

/***********************************************
Determine the imaginary potential beginning/end
************************************************/

  int MolPotenL = 0;
  int MolPotenR = 0;
  int found     = 0;
  while (found == 0)
  {
    if (abs(imag(MainCalc.Vgrid->element(IP.nx/2,MolPotenL))) == 0)
      found = 1;
    else
      MolPotenL++;
  }
  cout << "Left molecular absorbing potential ends at array index " << MolPotenL << "." << endl;

  MolPotenR = MolPotenL;
  while (found == 1)
  {
    if (abs(imag(MainCalc.Vgrid->element(IP.nx/2,MolPotenR))) != 0)
      found = 0;
    else
      MolPotenR++;
  }
  cout << "Right molecular absorbing potential begins at array index " << MolPotenR << "." << endl;



  return(0);
}
