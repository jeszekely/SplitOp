// #include <iostream>
// #include <cmath>
// #include <algorithm>
// #include <vector>
// #include <ctime>
#include "wvfxn.hpp"
#include "input_parser.hpp"
#include "splitop.hpp"
#include "junction.hpp"
#include "polynomial.hpp"
#include "chebyshev.hpp"

using namespace std;

int main(int argc, char const *argv[])
{
#if 1
  programInputs IP("inputs.json");

  SplitOp1D TestCalc(IP);
  transform(TestCalc.xgrid->data(),TestCalc.xgrid->data()+IP.nx,TestCalc.wvfxn->data(),[&](double x){return cplx(wvfxnElectron(x,IP));});
  TestCalc.wvfxn->normalize();
  TestCalc.initializeTDSE([&](double a){return cplx(0.0,-1.0)*absorbingPotential(a,IP)+Ve(a,IP);}, [&](double a){return HKinetic1D(a,IP.m_electron);});
  TestCalc.propagateStep(10);

  Chebyshev1D TestCheb(IP);
  transform(TestCheb.xgrid->data(),TestCheb.xgrid->data()+IP.nx,TestCheb.wvfxn->data(),[&](double x){return cplx(wvfxnElectron(x,IP));});
  TestCheb.wvfxn->normalize();
  TestCheb.dt *= 10;
  TestCheb.initializeTDSE([&](double a){return Ve(a,IP);}, [&](double a){return HKinetic1D(a,IP.m_electron);});
  TestCheb.propagateStep();
  printArrays(7000,cout,*TestCalc.xgrid,*TestCalc.Vgrid,*TestCalc.wvfxn,*TestCheb.wvfxn);
  //cout << TestCheb.wvfxn->getNorm() << endl << TestCalc.wvfxn->getNorm() << endl;

#endif

#if 0
  cout << "Hermite Test Case:" << endl;
  for (int ii = 0; ii < 10; ii++)
   cout << *Hermite<int>(ii);
  cout << endl << "Chebyshev Test Case:" << endl;
  for (int ii = 0; ii < 10; ii++)
  {
    cout << *Chebyshev<int>(ii);
    cout << Chebyshev<int>(ii)->eval(1.0) << endl << endl;
  }
 cout << endl;
  cout << *ClenshawChebyshevProp(4,1.0) << endl;
#endif

#if 0
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
  // for (int ii = 0; ii < IP.nx; ii++)
  //   for (int jj = 0; jj < IP.ny; jj++)
  //     MainCalc.wvfxn->element(ii,jj) = wvfxnElectron(MainCalc.xgrid->element(ii),IP)*MoleStates(0,jj);
  // MainCalc.wvfxn->normalize();
  // MainCalc.propagateStep(10);

#endif

#if 0
  ElecCalc.initializeTDSE([&](double a){return ElecParab(a,IP.requil,IP);}, [&](double a){return HKinetic1D(a,IP.m_electron);});

  Array1D<double> TestCheb(ElecCalc.PotenOp->Nx(),0.0,0.0);
  Array1D<cplx> TestCheb2(ElecCalc.PotenOp->Nx(),0.0,0.0);

//get min/max for chebychev calc
  for (int ii = 0; ii < ElecCalc.Vgrid->Nx(); ++ii)
    TestCheb(ii) = real(ElecCalc.Vgrid->element(ii));
  for (int ii = 0; ii < ElecCalc.Vgrid->Nx(); ++ii)
    TestCheb2(ii) = ElecCalc.Vgrid->element(ii);

  auto result = minmax_element(&TestCheb(0),&TestCheb(0)+TestCheb.Nx());
  double min  = *result.first;
  double max  = *result.second;
  cout << "Min: " << min << endl;
  cout << "Max: " << max << endl;

  double alpha = 0.5*IP.dt*(max - min);
  double dE    = 0.5*(max + min);

  for (int ii = 0; ii < TestCheb.Nx(); ++ii)
    TestCheb2(ii) -= dE; //Need an identity if this is a matrix
  TestCheb2.scale(2.0/(max-min));

  shared_ptr<polynomial<cplx>> ChebCoeffs = ClenshawChebyshevProp(10,alpha);
  cplx coeff = exp(cplx(0.0,-1.0)*(min+alpha));

  Array1D<cplx> WorkArray(ElecCalc.PotenOp->Nx(),0.0,0.0);
  cout << "Using " << ChebCoeffs->vals->size() << " terms in the expansion." << endl;

  WorkArray = ChebCoeffs->eval(TestCheb2);
  WorkArray.scale(coeff);

   for (int ii = 0; ii < ElecCalc.PotenOp->size(); ii++)
     cout << ii << " " << ElecCalc.PotenOp->element(ii) << " " << WorkArray(ii) << endl;

#endif

  return(0);
}
