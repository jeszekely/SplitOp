#include "surrogate.hpp"
#include "chebyshev.hpp"
#include "utilities.hpp"

using namespace std;

wvfxn1DArray::wvfxn1DArray(const int nn, const double xstep, programInputs &IP) : Array2D<cplx>(IP.nx,nn,xstep,1.0)
{
  mass = IP.m_electron;
  hbar = IP.hbar;
}

double wvfxn1DArray::getNorm(const int nn)
{
  wvfxn1D S(nx,0.0,xstep);
  copy_n(&element(0,nn),nx,S.data());
  auto T = S | S;

  return real(T.integrate_rect());
}

void wvfxn1DArray::normalize(const int nn)
{
  double norm = getNorm(nn);
  assert (norm > 0);
  double fact = 1.0/sqrt(norm);
  std::for_each(&element(0,nn),&element(0,nn)+nx, [&fact](cplx &p){p*=fact;});
}

double wvfxn1DArray::flux(const int nn, const int xx)
{
  cplx psi_star = conj(element(nn,xx));
  cplx psi_deriv = cplx(0,-1.0)*deriv_16_x(nn,xx);
  return real(psi_deriv*psi_star*hbar/mass);
}

double wvfxn1DArray::m() {return mass;}
double wvfxn1DArray::hb() {return hbar;}


SplitOp1DArray::SplitOp1DArray(string filename) : Ins(programInputs(filename))
{
//Define constants
  nx       = Ins.nx;
  n        = Ins.nphonon;
  nthreads = Ins.procs;
  xmin     = Ins.xmin;
  xmax     = Ins.xmax;
  xstep    = (xmax-xmin)/nx;
  pstep    = 2.0*M_PI/(nx*xstep);
  dt       = Ins.dt;
  simtime  = 0.0;
  runtime  = Ins.runtime;

//Define Arrays
  xgrid = make_shared<Array1D<double>>(nx,xmin,xstep);
  xgrid->fill_array();

  pgrid = make_shared<Array1D<double>>(nx,0.0,pstep);
  pgrid->fill_array();
  for (int ii = nx/2; ii < nx; ii++)
    pgrid->element(ii) = -1.0*pgrid->element(nx-ii);

  Tgrid   = make_shared<Array1D<double>>(nx,xmin,pstep);
  KinetOp = make_shared<Array1D<cplx>>(nx,xmin,pstep);
  wvfxn   = make_shared<wvfxn1DArray>(n, xstep, Ins);

  fftw_init_threads(); //initialize SMP
  fftw_plan_with_nthreads(nthreads); //use #(procs) processors

  for (int ii = 0; ii < n; ii++)
  {
    forplanArray.push_back(fftw_plan_dft_1d(nx,reinterpret_cast<fftw_complex*>(&wvfxn->element(0,ii)),reinterpret_cast<fftw_complex*>(&wvfxn->element(0,ii)),FFTW_FORWARD,FFTW_MEASURE));
    backplanArray.push_back(fftw_plan_dft_1d(nx,reinterpret_cast<fftw_complex*>(&wvfxn->element(0,ii)),reinterpret_cast<fftw_complex*>(&wvfxn->element(0,ii)),FFTW_BACKWARD,FFTW_MEASURE));
  }
}

SplitOp1DArray::~SplitOp1DArray()
{
  for (fftw_plan &plan : forplanArray)
    fftw_destroy_plan(plan);
  for (fftw_plan &plan : backplanArray)
    fftw_destroy_plan(plan);
}

void SplitOp1DArray::initializeSurrogate()
{
  transform(pgrid->data(),pgrid->data()+nx,Tgrid->data(),[&](double a){return HKinetic1D(a,Ins.m_electron);});
  transform(Tgrid->data(),Tgrid->data()+nx,KinetOp->data(),[&](double a){return exp(cplx(0.0,-0.5)*a*Ins.dt);});

  for (int ii = 0; ii < nx; ii++)
    HInteractions.push_back(make_shared<matrixComp>(n,n));

  //get the phonon eigenstates
  programInputs MoleIP = Ins;
  MoleIP.nx            = Ins.ny;
  MoleIP.xmin          = Ins.ymin;
  MoleIP.xmax          = Ins.ymax;
  MoleIP.m_electron    = Ins.m_C60;
  SplitOp1D MoleCalc(MoleIP);
  MoleCalc.initializeITP([&](double a){return VM(a,Ins);}, [&](double a){return HKinetic1D(a,Ins.m_C60);});
  Array2D <cplx> MoleStates(Ins.ny,Ins.nphonon,Ins.ymin,MoleCalc.xstep);
  MoleCalc.propagateITP(MoleStates);
  vector<double> PhononEnergies;
  for (int ii = 0; ii < n; ii++)
  {
    copy_n(&MoleStates(0,ii),Ins.ny,MoleCalc.wvfxn->data());
    cout << ii << " " <<  MoleCalc.getEnergy() << endl;
    PhononEnergies.push_back(MoleCalc.getEnergy());
  }
  if (Ins.phononstates == 1)
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

  vector<shared_ptr<wvfxn1D>> PhononStates;
  for (int ii = 0; ii < n; ii++)
    PhononStates.push_back(make_shared<wvfxn1D>(MoleCalc.xgrid->size(),MoleCalc.xgrid->init(),MoleCalc.xgrid->step()));
  for (int ii = 0; ii < n; ii++)
    copy_n(&MoleStates(0,ii),MoleCalc.xgrid->size(),PhononStates[ii]->data());

  //Fill the interaction matrices with subspace coupling terms
  Array1D <cplx> HInt(MoleCalc.xgrid->size(),MoleCalc.xgrid->init(),MoleCalc.xgrid->step());
  for (int ii = 0; ii < nx; ii++)
  {
    for (int jj = 0; jj < MoleCalc.xgrid->size(); jj++)
      HInt(jj) = VeM(xgrid->element(ii), MoleCalc.xgrid->element(jj), Ins);
    Array1D <cplx> Scratch(MoleCalc.xgrid->size(),MoleCalc.xgrid->init(),MoleCalc.xgrid->step());
    for (int jj = 0; jj < n; jj++)
    {
      for (int kk = 0; kk < jj; kk++)
      {
        Scratch = *PhononStates[kk] | *PhononStates[jj];
        Scratch *= HInt;
        HInteractions[ii]->element(jj,kk) = Scratch.integrate_rect();
        HInteractions[ii]->element(kk,jj) = HInteractions[ii]->element(jj,kk);
      }
    }
  }
  //print some of the coupling matrix elements
  ofstream CoupingFile;
    CoupingFile.open("output_data/CoupingResults.txt");
    for (int kk = 0; kk < xgrid->size(); kk++)
      CoupingFile << xgrid->element(kk) << " " << real(HInteractions[kk]->element(0,1)) << " " << real(HInteractions[kk]->element(1,2)) << " " << real(HInteractions[kk]->element(0,2)) << endl;
    CoupingFile.close();

  //Add other energy terms to interaction operators
  for (int ii = 0; ii < xgrid->size(); ii++)
  {
    for (int jj = 0; jj < n; jj++)
      HInteractions[ii]->element(jj,jj) += ( PhononEnergies[jj] + Ve(xgrid->element(ii),Ins) );
  }
  //Exponentiate the interaction matrices
  for (int ii = 0; ii < xgrid->size(); ii++)
    exponentiateCheb(*HInteractions[ii],dt);
}

void SplitOp1DArray::propagateStep(int nn)
{
  //Apply kinetic operators
  for (fftw_plan &plan : forplanArray)
    fftw_execute(plan);
  for (int ii = 0; ii<n; ii++)
    transform(&wvfxn->element(0,ii),&wvfxn->element(0,ii)+nx,KinetOp->data(),&wvfxn->element(0,ii),multiplies<cplx>());
  for (fftw_plan &plan : backplanArray)
    fftw_execute(plan);
  wvfxn->scale(1.0/double(nx));

  //Apply potential/interaction operators
  for (int ii = 0; ii < xgrid->size(); ii++)
  {
    matrixComp scratch(n,1);
    for (int jj = 0; jj < n; jj++)
      scratch(jj,0) = wvfxn->element(ii,jj);
    zgemv_("N",n,n,cplx(1.0,0.0),HInteractions[ii]->data(),n,scratch.data(),1,cplx(0.0,0.0),&wvfxn->element(ii,0),xgrid->size());
  }

  //Apply kinetic operators
  for (fftw_plan &plan : forplanArray)
    fftw_execute(plan);
  for (int ii = 0; ii<n; ii++)
    transform(&wvfxn->element(0,ii),&wvfxn->element(0,ii)+nx,KinetOp->data(),&wvfxn->element(0,ii),multiplies<cplx>());
  for (fftw_plan &plan : backplanArray)
    fftw_execute(plan);

  wvfxn->scale(1.0/double(nx));

  simtime += nn*dt;
}

void SplitOp1DArray::printWavefunction(std::string filename)
{
  ofstream out;
  out.open(filename);
  for (int ii = 0; ii < nx; ii++)
  {
    out << xgrid->element(ii) << " ";
    for (int jj = 0; jj < n; jj++)
      out << abs(wvfxn->element(ii,jj)) << " ";
    out << endl;
  }
  out.close();
}

