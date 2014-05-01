#include "chebyshev.hpp"
#include "utilities.hpp"

using namespace std;

double ChebyshevCoeff(int nn, double alpha)
{
  //cout << nn << " " << alpha << " " << boost::math::cyl_bessel_j(nn,alpha) << endl;
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
  polyterms= 200;

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
  //result = minmax_element(&Tgrid->element(0),&Tgrid->element(0)+Tgrid->Nx());
  Tmin  = 0.0;//*result.first;
  Tmax  = pow(M_PI/xstep,2)*0.5;//*result.second;
}

void Chebyshev1D::ApplyHNorm(shared_ptr<wvfxn1D> W)
{
  //cout << W->getNorm()<< endl;
  //Store a copy of the wavefunction for later
  wvfxn1D o(*W);
  auto W1 = make_shared<wvfxn1D>(*W);
  fftw_plan forplan  = fftw_plan_dft_1d(nx,reinterpret_cast<fftw_complex*>(W1->data()),reinterpret_cast<fftw_complex*>(W1->data()),FFTW_FORWARD,FFTW_MEASURE);
  fftw_plan backplan = fftw_plan_dft_1d(nx,reinterpret_cast<fftw_complex*>(W1->data()),reinterpret_cast<fftw_complex*>(W1->data()),FFTW_BACKWARD,FFTW_MEASURE);
  //copy data back to W in case fftw plans modified the data
  copy_n(o.data(),nx,W1->data());
  double dE = Vmax+Tmax-Vmin-Tmin;
  //o.scale(-1.0*(1.0+2.0*(Vmin+Tmin)/dE));

  //Apply normalized potential operator
  transform(W->data(),W->data()+nx,Vgrid->data(),W->data(),multiplies<cplx>());
  for (int ii = 0; ii < W->size(); ii++)
    W->element(ii) -= Vmin;
  W->scale(2.0/dE);

  //Apply normalized kinetic operator
  fftw_execute(forplan);
  transform(W1->data(),W1->data()+nx,Tgrid->data(),W1->data(),multiplies<cplx>());
  fftw_execute(backplan);
  W1->scale(1.0/double(nx));
  W1->scale(2.0/dE);

  //transform(W->data(),W->data()+nx,W1->data(),W->data(),plus<cplx>());
  zaxpy_(W->size(), cplx(1.0,0.0), W1->data(), 1, W->data(), 1);
  //*W += o;
  zaxpy_(W->size(), cplx(-1.0,0.0), o.data(), 1, W->data(), 1);
  fftw_destroy_plan(forplan);
  fftw_destroy_plan(backplan);

  //Combine Results
  //cout << W->getNorm() << endl << endl;
}

void Chebyshev1D::propagateStep()
{
  double alpha = 0.5*dt*(Vmax+Tmax-Vmin-Tmin);
  shared_ptr<wvfxn1D> phi1 = make_shared<wvfxn1D>(*wvfxn);
  ApplyHNorm(phi1);
  phi1->scale(cplx(0.0,-1.0));
  //shared_ptr<wvfxn1D> bkp1 = make_shared<wvfxn1D>(nx,0.0,0.0);
  //shared_ptr<wvfxn1D> bkp2 = make_shared<wvfxn1D>(nx,0.0,0.0);
  //shared_ptr<wvfxn1D> bk   = make_shared<wvfxn1D>(nx,0.0,0.0);
  //int kk                   = polyterms;
  shared_ptr<wvfxn1D> phinm1 = make_shared<wvfxn1D>(*phi1);
  shared_ptr<wvfxn1D> phinm2 = make_shared<wvfxn1D>(*wvfxn);
  shared_ptr<wvfxn1D> phin   = make_shared<wvfxn1D>(*wvfxn);
  phin->zero();
  double ak;
  //while (true)
  //Apply first two terms of expansion
  wvfxn->scale(ChebyshevCoeff(0,alpha));
  ak = ChebyshevCoeff(1,alpha);

  zaxpy_(wvfxn->size(), ak, phi1->data(), 1, wvfxn->data(), 1);

  for (int kk = 2; kk <= polyterms; kk++)
  {
    ak = ChebyshevCoeff(kk,alpha);
    copy_n(phinm2->data(),nx,phin->data()); //don't need data in phinm2 after this step of the loop
    copy_n(phinm1->data(),nx,phinm2->data());
    //*phinm2 = *phinm1;
    ApplyHNorm(phinm1);
    phinm1->scale(cplx(0.0,-2.0));
    zaxpy_(phin->size(), cplx(1.0,0.0), phinm1->data(), 1, phin->data(), 1);
    copy_n(phin->data(),nx,phinm1->data());
    //*phinm1 = *phin;
    zaxpy_(wvfxn->size(), ak, phin->data(), 1, wvfxn->data(), 1);
    phin->zero();
    // cout << kk << " " << wvfxn->getNorm() << " " << ak << " " << ak*ak*phinm1->getNorm() << endl;

    // ak = ChebyshevCoeff(kk,alpha);
    // phin = make_shared<wvfxn1D>(*phinm2);
    // phinm2 = make_shared<wvfxn1D>(*phinm1);
    // ApplyHNorm(phinm2);
    // phinm2->scale(cplx(0.0,-2.0));
    // transform(phin->data(),phin->data()+nx,phinm2->data(),phin->data(),plus<cplx>());
    // transform(wvfxn->data(),wvfxn->data()+nx,phin->data(),wvfxn->data(),[&](cplx &a, cplx &b){return a+b*ak;});
    // *phinm2 = *phinm1;
    // *phinm1 = *phin;
        //cout << *phinm2 << *phinm1 << *phin << endl;

    // ApplyHNorm(bkp1);
    // bkp1->scale(cplx(0.0,-2.0));
    // transform(bkp1->data(),bkp1->data()+nx,bkp2->data(),bk->data(),plus<cplx>());
    // for_each(bk->data(),bk->data()+nx,[&](cplx &n){n+=ak;});
    // if (kk == 0) break;
    // *bkp2 = *bkp1;
    // *bkp1 = *bk;
    // --kk;
  }
  // auto wvfxn1 = make_shared<wvfxn1D>(*phi1);
  // transform(wvfxn1->data(),wvfxn1->data()+nx,bkp1->data(),wvfxn1->data(),multiplies<cplx>());
  // auto wvfxn2 = make_shared<wvfxn1D>(phi0);
  // transform(wvfxn2->data(),wvfxn2->data()+nx,bkp2->data(),wvfxn2->data(),multiplies<cplx>());
  // auto wvfxn3 = make_shared<wvfxn1D>(phi0);
  // for_each(wvfxn3->data(),wvfxn3->data()+nx,[&](cplx &n){n*=ak;});
  // wvfxn->zero();
  // *wvfxn+= *wvfxn1;// = make_shared<wvfxn1D>(*wvfxn1 + *wvfxn2 + *wvfxn3);
  // *wvfxn+= *wvfxn2;
  // *wvfxn+= *wvfxn3;
  cplx coeff = exp(cplx(0.0,-1.0)*(Vmin+Tmin+alpha));
  wvfxn->scale(coeff);

  // polynomial<cplx> phi0(1);
  // phi0(0)  = 1;
  // polynomial<cplx> phi1(2);
  // phi1(1)  = cplx(0.0,-1.0);
  // polynomial<cplx> bkp1(1);
  // polynomial<cplx> bkp2(1);
  // polynomial<cplx> alpha(2);
  //alpha(1) = cplx(0.0,-2.0);
  //polynomial<cplx> beta(1);
  //beta(0)  = 1; //-1 for normal chebyshev iterative definition
  // int kk   = nn;
  // polynomial<cplx> bk(1);
  // polynomial<cplx> ak(1);
  // polynomial<cplx> a0(1);
  // while (true)
  // {
  //   ak(0) = ChebyshevCoeff(kk,a);
  //   bk    = alpha*bkp1 + beta*bkp2 + ak;
  //   //cout << kk << ":" << endl << bk << bkp1 <<  bkp2 << phi1*bkp1+beta*phi0*bkp2 << endl << endl;
  //   if (kk == 0) break;
  //   bkp2  = bkp1;
  //   bkp1  = bk;
  //   --kk;
  // }
  // shared_ptr<polynomial<cplx>> S = make_shared<polynomial<cplx>>(phi1*bkp1+beta*phi0*bkp2+phi0*ak);
  // return S;

}