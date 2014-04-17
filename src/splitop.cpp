#include "splitop.hpp"

using namespace std;

SplitOp1D::SplitOp1D(programInputs &IP)
{
//Define constants
  nx = IP.nx;
  nthreads = IP.procs;
  xmin = IP.xmin;
  xmax = IP.xmax;
  xstep = (xmax-xmin)/nx;
  pstep = 2.0*M_PI/(nx*xstep);
  dt = IP.dt;
  simtime = 0.0;
  runtime = IP.runtime;

//Define Arrays
  xgrid = make_shared<Array1D<double>>(nx,xmin,xstep);
  xgrid->fill_array();

  pgrid = make_shared<Array1D<double>>(nx,0.0,pstep);
  pgrid->fill_array();
  for (int ii = nx/2; ii < nx; ii++)
    pgrid->element(ii) = -1.0*pgrid->element(nx-ii);

 	Vgrid   = make_shared<Array1D<cplx>>(nx,xmin,xstep);
  Tgrid   = make_shared<Array1D<double>>(nx,xmin,xstep);
  KinetOp = make_shared<Array1D<cplx>>(nx,xmin,xstep);
  PotenOp = make_shared<Array1D<cplx>>(nx,xmin,xstep);
  wvfxn   = make_shared<wvfxn1D>(nx,xmin,xstep);

  fftw_init_threads(); //initialize SMP
  fftw_plan_with_nthreads(nthreads); //use #(procs) processors

  forplan  = fftw_plan_dft_1d(nx,(fftw_complex *)wvfxn->data(),(fftw_complex *)wvfxn->data(),FFTW_FORWARD,FFTW_MEASURE);
  backplan = fftw_plan_dft_1d(nx,(fftw_complex *)wvfxn->data(),(fftw_complex *)wvfxn->data(),FFTW_FORWARD,FFTW_MEASURE);
}

void SplitOp1D::initializeTDSE(std::function<cplx(double)> fV, std::function<double(double)> fT)
{
	transform(xgrid->data(),xgrid->data()+nx,Vgrid->data(),fV);
	transform(pgrid->data(),pgrid->data()+nx,Tgrid->data(),fT);
	transform(Vgrid->data(),Vgrid->data()+nx,PotenOp->data(),[&](cplx a){return exp(cplx(0.0,-1.0)*a*dt);});
}

