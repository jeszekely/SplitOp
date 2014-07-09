#include <ctime>
#include "wvfxn.hpp"
#include "input_parser.hpp"
#include "splitop.hpp"
#include "junction.hpp"
#include "polynomial.hpp"
#include "chebyshev.hpp"
#include "surrogate.hpp"

using namespace std;

int main(int argc, char const *argv[])
{
  programInputs IP("inputs.json");

#if 1
  SplitOp1DArray TestSurr("inputs.json");
  TestSurr.initializeSurrogate();
  transform(TestSurr.xgrid->data(),TestSurr.xgrid->data()+IP.nx,&TestSurr.wvfxn->element(0,0),[&](double x){return cplx(wvfxnElectron(x,IP));});
  TestSurr.wvfxn->normalize(0);
  for (int ii = 0; ii < 800; ii++)
  {
    cout << ii << " " << TestSurr.wvfxn->getNorm(0) + TestSurr.wvfxn->getNorm(1) + TestSurr.wvfxn->getNorm(2) << endl;
    TestSurr.propagateStep();
  }
  TestSurr.printWavefunction("output_data/FinalWavefunction.txt");


#endif

#if 0
  clock_t t1, t2;
  SplitOp1D TestCalc(IP);
  transform(TestCalc.xgrid->data(),TestCalc.xgrid->data()+IP.nx,TestCalc.wvfxn->data(),[&](double x){return cplx(wvfxnElectron(x,IP));});
  TestCalc.wvfxn->normalize();
  TestCalc.initializeTDSE([&](double a){return cplx(0.0,-1.0)*absorbingPotential(a,IP)+Ve(a,IP);}, [&](double a){return HKinetic1D(a,IP.m_electron);});
  t1 = clock();
  TestCalc.propagateStep(100);
  t2 = clock();
  double splitopTime = (double(t2) - double(t1))/CLOCKS_PER_SEC;
  cout << "Split Operator Time:        " << setw(20) << setprecision(8) << splitopTime << endl;

  Chebyshev1D TestCheb(IP);
  transform(TestCheb.xgrid->data(),TestCheb.xgrid->data()+IP.nx,TestCheb.wvfxn->data(),[&](double x){return cplx(wvfxnElectron(x,IP));});
  TestCheb.wvfxn->normalize();
  TestCheb.dt = 10;
  TestCheb.initializeTDSE([&](double a){return Ve(a,IP);}, [&](double a){return HKinetic1D(a,IP.m_electron);});
  t1 = clock();
  TestCheb.propagateStep();
  t2 = clock();
  double chebyTime = (double(t2) - double(t1))/CLOCKS_PER_SEC;
  cout << "Chebyshev Propagation Time: " << setw(20) << setprecision(8) << chebyTime << endl;

  //printArrays(7000,cout,*TestCalc.xgrid,*TestCalc.Vgrid,*TestCalc.wvfxn,*TestCheb.wvfxn);
  //cout << TestCheb.wvfxn->getNorm() << endl << TestCalc.wvfxn->getNorm() << endl;

#endif

#if 0

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
Determine the imaginary potentials begin/end
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

  int ElecPotenL = 0;
  int ElecPotenR = 0;
  found          = 0;
  while (found == 0)
  {
    if (abs(imag(MainCalc.Vgrid->element(ElecPotenL,IP.ny/2))) == 0)
      found = 1;
    else
      ElecPotenL++;
  }
  cout << "Left electronic absorbing potential ends at array index " << ElecPotenL << "." << endl;

  ElecPotenR = ElecPotenL;
  while (found == 1)
  {
    if (abs(imag(MainCalc.Vgrid->element(ElecPotenR,IP.ny/2))) != 0)
      found = 0;
    else
      ElecPotenR++;
  }
  cout << "Right electronic absorbing potential begins at array index " << ElecPotenR << "." << endl;


/***********************************************
Determine the junction bounds
************************************************/

  int zRjunc = 0;
  int zLjunc = 0;
  for (int ii = 0; ii < IP.nx; ii++)
  {
    if (MainCalc.xgrid->element(ii) >= IP.zl)
    {
      zLjunc = ii;
      break;
    }
  }
  for (int ii = 0; ii < IP.nx; ii++)
  {
    if (MainCalc.xgrid->element(ii) >= IP.zr)
    {
      zRjunc = ii;
      break;
    }
  }
  cout << "The junction spans from index " << zLjunc << " to index " << zRjunc << "." << endl;

/***********************************************
Define initial 2D wavefunction
************************************************/
  for (int ii = 0; ii < IP.nx; ii++)
    for (int jj = 0; jj < IP.ny; jj++)
      MainCalc.wvfxn->element(ii,jj) = wvfxnElectron(MainCalc.xgrid->element(ii),IP)*MoleStates(0,jj);
  MainCalc.wvfxn->normalize();
  MainCalc.propagateStep(10);

#endif

  return(0);
}
