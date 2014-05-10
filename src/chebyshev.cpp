#include "chebyshev.hpp"
#include "utilities.hpp"

using namespace std;

double ChebyshevCoeff(int nn, double alpha)
{
  if (nn == 0) return boost::math::cyl_bessel_j(nn,alpha);
  else return 2.0*boost::math::cyl_bessel_j(nn,alpha);
}

Chebyshev1D::Chebyshev1D(programInputs &IP)
{
//Define constants
  nx       = IP.nx;
  nthreads = IP.procs;
  xmin     = IP.xmin;
  xmax     = IP.xmax;
  xstep    = (xmax-xmin)/nx;
  pstep    = 2.0*M_PI/(nx*xstep);
  dt       = IP.dt;
  itpdt    = 1.0;
  simtime  = 0.0;
  runtime  = IP.runtime;

//Define Arrays
  xgrid = make_shared<Array1D<double>>(nx,xmin,xstep);
  xgrid->fill_array();

  pgrid = make_shared<Array1D<double>>(nx,0.0,pstep);
  pgrid->fill_array();
  for (int ii = nx/2; ii < nx; ii++)
    pgrid->element(ii) = -1.0*pgrid->element(nx-ii);

  Vgrid   = make_shared<Array1D<double>>(nx,xmin,xstep);
  Tgrid   = make_shared<Array1D<double>>(nx,xmin,pstep);
  wvfxn   = make_shared<wvfxn1D>(nx,xmin,xstep,IP.hbar,IP.m_electron);

  fftw_init_threads(); //initialize SMP
  fftw_plan_with_nthreads(nthreads); //use #(procs) processors
}

void Chebyshev1D::initializeTDSE(std::function<double(double)> fV, std::function<double(double)> fT)
{
  transform(xgrid->data(),xgrid->data()+nx,Vgrid->data(),fV);
  transform(pgrid->data(),pgrid->data()+nx,Tgrid->data(),fT);
  auto result = minmax_element(&Vgrid->element(0),&Vgrid->element(0)+Vgrid->Nx());
  Vmin  = *result.first;
  Vmax  = *result.second;
  Tmin  = 0.0;//*result.first;
  Tmax  = pow(M_PI/xstep,2)*0.5;//*result.second;
}

void Chebyshev1D::ApplyHNorm(shared_ptr<wvfxn1D> W)
{
  //Store a copy of the wavefunction for later
  wvfxn1D o(*W);
  auto W1            = make_shared<wvfxn1D>(*W);
  fftw_plan forplan  = fftw_plan_dft_1d(nx,reinterpret_cast<fftw_complex*>(W1->data()),reinterpret_cast<fftw_complex*>(W1->data()),FFTW_FORWARD,FFTW_MEASURE);
  fftw_plan backplan = fftw_plan_dft_1d(nx,reinterpret_cast<fftw_complex*>(W1->data()),reinterpret_cast<fftw_complex*>(W1->data()),FFTW_BACKWARD,FFTW_MEASURE);
  //copy data back to W in case fftw plans modify the data
  copy_n(o.data(),nx,W1->data());
  double dE          = Vmax+Tmax-Vmin-Tmin;

  //Apply normalized potential operator
  transform(W->data(),W->data()+nx,Vgrid->data(),W->data(),multiplies<cplx>());
  zaxpy_(W->size(), cplx(-1.0*Vmin), o.data(), 1, W->data(), 1);

  //Apply normalized kinetic operator
  fftw_execute(forplan);
  transform(W1->data(),W1->data()+nx,Tgrid->data(),W1->data(),multiplies<cplx>());
  fftw_execute(backplan);
  W1->scale(1.0/double(nx));

  zaxpy_(W->size(), cplx(1.0,0.0), W1->data(), 1, W->data(), 1);
  W->scale(2.0/dE);
  zaxpy_(W->size(), cplx(-1.0,0.0), o.data(), 1, W->data(), 1);
  fftw_destroy_plan(forplan);
  fftw_destroy_plan(backplan);
}

void Chebyshev1D::propagateStep()
{
  double alpha               = 0.5*dt*(Vmax+Tmax-Vmin-Tmin);
  shared_ptr<wvfxn1D> phi1   = make_shared<wvfxn1D>(*wvfxn);
  ApplyHNorm(phi1);
  phi1->scale(cplx(0.0,-1.0));
  shared_ptr<wvfxn1D> phinm1 = make_shared<wvfxn1D>(*phi1);
  shared_ptr<wvfxn1D> phinm2 = make_shared<wvfxn1D>(*wvfxn);
  shared_ptr<wvfxn1D> phin   = make_shared<wvfxn1D>(*wvfxn);
  phin->zero();
  double ak;

  //Apply first two terms of expansion
  wvfxn->scale(ChebyshevCoeff(0,alpha));

  ak = ChebyshevCoeff(1,alpha);
  zaxpy_(wvfxn->size(), ak, phi1->data(), 1, wvfxn->data(), 1);
  polyterms = 200+int((Vmax+Tmax-Vmin-Tmin)*dt*0.5);
  for (int kk = 2; kk <= polyterms; kk++)
  {
    ak = ChebyshevCoeff(kk,alpha);
    copy_n(phinm2->data(),nx,phin->data()); //don't need data in phinm2 after this step of the loop
    copy_n(phinm1->data(),nx,phinm2->data());

    ApplyHNorm(phinm1);
    zaxpy_(phin->size(), cplx(0.0,-2.0), phinm1->data(), 1, phin->data(), 1);
    copy_n(phin->data(),nx,phinm1->data());

    zaxpy_(wvfxn->size(), ak, phin->data(), 1, wvfxn->data(), 1);
    phin->zero();
    //cout << kk << " " << wvfxn->getNorm() << " " << ak << " " << ak*ak*phinm1->getNorm() << endl;
  }
  //cout << polyterms << endl;
  cplx coeff = exp(cplx(0.0,-1.0)*(Vmin+Tmin+alpha));
  wvfxn->scale(coeff);
}

void exponentiateCheb(matrixComp &H, double dt)
{
  vector<double> HEvals(H.nr(),0.0);
  H.getEigvals(HEvals.data());
  double min = HEvals[0];
  double max = HEvals[H.nr()-1];

  //Form normalized H
  double dE = max - min;
  matrixComp Hn(H);
  Hn.makeIdentity();
  Hn.scale(min-0.5*dE);
  H         -= Hn;
  H.scale(2.0/dE);

//Set up matricies for expansion
  double alpha = 0.5*dt*dE;
  matrixComp phi1(Hn);
  phi1.scale(cplx(0.0,-1.0));
  Hn           = H;

  matrixComp phinm1(phi1);
  matrixComp phinm2(H);
  matrixComp phin(H);
  phin.zero();
  double ak;

  H.makeIdentity();
  H.scale(ChebyshevCoeff(0,alpha));
  ak            = ChebyshevCoeff(1,alpha);
  H             += phi1*ak;
  int polyterms = 10+int(dE*dt*0.5);

  for (int kk = 2; kk <= polyterms; kk++)
  {
    ak     = ChebyshevCoeff(kk,alpha);
    phin   = phinm2;
    phinm2 = phinm1;

    phinm1 *= Hn;
    phinm1.scale(cplx(0.0,-2.0));
    phin   += phinm1;
    phinm1 = phin;
    H      += phin*ak;
    phin.zero();
  }
  H.scale(exp(cplx(0.0,-1.0)*(min+alpha)));
}


