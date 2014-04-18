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
  itpdt = 1.0;
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
  Tgrid   = make_shared<Array1D<double>>(nx,xmin,pstep);
  KinetOp = make_shared<Array1D<cplx>>(nx,xmin,pstep);
  PotenOp = make_shared<Array1D<cplx>>(nx,xmin,xstep);
  wvfxn   = make_shared<wvfxn1D>(nx,xmin,xstep);

  fftw_init_threads(); //initialize SMP
  fftw_plan_with_nthreads(nthreads); //use #(procs) processors

  forplan  = fftw_plan_dft_1d(nx,reinterpret_cast<fftw_complex*>(wvfxn->data()),reinterpret_cast<fftw_complex*>(wvfxn->data()),FFTW_FORWARD,FFTW_MEASURE);
  backplan = fftw_plan_dft_1d(nx,reinterpret_cast<fftw_complex*>(wvfxn->data()),reinterpret_cast<fftw_complex*>(wvfxn->data()),FFTW_BACKWARD,FFTW_MEASURE);
}

SplitOp1D::~SplitOp1D(){fftw_destroy_plan(forplan); fftw_destroy_plan(backplan);}

void SplitOp1D::initializeTDSE(std::function<cplx(double)> fV, std::function<double(double)> fT)
{
	transform(xgrid->data(),xgrid->data()+nx,Vgrid->data(),fV);
	transform(pgrid->data(),pgrid->data()+nx,Tgrid->data(),fT);
	transform(Vgrid->data(),Vgrid->data()+nx,PotenOp->data(),[&](cplx a){return exp(cplx(0.0,-1.0)*a*dt);});
	transform(Tgrid->data(),Tgrid->data()+nx,KinetOp->data(),[&](double a){return exp(cplx(0.0,-0.5)*a*dt);});
}

void SplitOp1D::propagateStep()
{
	//Apply kinetic operator
	fftw_execute(forplan);
	transform(wvfxn->data(),wvfxn->data()+nx,KinetOp->data(),wvfxn->data(),multiplies<cplx>());
	fftw_execute(backplan);
	wvfxn->scale(1.0/double(nx));

	//Apply potential operator
	transform(wvfxn->data(),wvfxn->data()+nx,PotenOp->data(),wvfxn->data(),multiplies<cplx>());

	//Apply kinetic operator
	fftw_execute(forplan);
	transform(wvfxn->data(),wvfxn->data()+nx,KinetOp->data(),wvfxn->data(),multiplies<cplx>());
	fftw_execute(backplan);
	wvfxn->scale(1.0/double(nx));

	simtime += dt;
}

void SplitOp1D::propagateNSteps(int nn)
{
	fftw_execute(forplan);
	transform(wvfxn->data(),wvfxn->data()+nx,KinetOp->data(),wvfxn->data(),multiplies<cplx>());
	fftw_execute(backplan);
	wvfxn->scale(1.0/nx);

	transform(wvfxn->data(),wvfxn->data()+nx,PotenOp->data(),wvfxn->data(),multiplies<cplx>());

	for (int ii = 0; ii < nn-1; ii++)
	{
		fftw_execute(forplan);
		transform(wvfxn->data(),wvfxn->data()+nx,KinetOp->data(),wvfxn->data(),[](cplx a, cplx b){return a*pow(b,2);});
		fftw_execute(backplan);
		wvfxn->scale(1.0/nx);
		transform(wvfxn->data(),wvfxn->data()+nx,PotenOp->data(),wvfxn->data(),multiplies<cplx>());
	}

	fftw_execute(forplan);
	transform(wvfxn->data(),wvfxn->data()+nx,KinetOp->data(),wvfxn->data(),multiplies<cplx>());
	fftw_execute(backplan);
	wvfxn->scale(1.0/nx);
	simtime += nn*dt;
}

void SplitOp1D::initializeITP(std::function<cplx(double)> fV, std::function<double(double)> fT)
{
	transform(xgrid->data(),xgrid->data()+nx,Vgrid->data(),fV);
	transform(pgrid->data(),pgrid->data()+nx,Tgrid->data(),fT);
	transform(Vgrid->data(),Vgrid->data()+nx,PotenOp->data(),[&](cplx a){return exp(-1.0*a*itpdt);});
	transform(Tgrid->data(),Tgrid->data()+nx,KinetOp->data(),[&](double a){return exp(-0.5*a*itpdt);});
}

void SplitOp1D::updateITP()
{
	transform(Vgrid->data(),Vgrid->data()+nx,PotenOp->data(),[&](cplx a){return exp(-1.0*a*itpdt);});
	transform(Tgrid->data(),Tgrid->data()+nx,KinetOp->data(),[&](double a){return exp(-0.5*a*itpdt);});
}

void SplitOp1D::propagateITP(Array2D<cplx> &States, double tolerance)
{
	assert(States.Ny() >= 1 && States.Nx() == nx);
	double conv_prev 	= 1.0;
	double conv_curr 	= 1.0;
	double conv_ratio 	= conv_curr/conv_prev;

	int nv = States.Ny();
	int steps = 0;
	vector<shared_ptr<wvfxn1D>> TestVecs;
	for (int ii = 0; ii < nv; ii++)
		TestVecs[ii] = make_shared<wvfxn1D>(nx,0.0,xstep);
	initializeGuess(TestVecs);

//Get the ground state vector
	copy_n(TestVecs[0]->data(),nx,wvfxn->data());
	wvfxn->normalize();
	while (conv_ratio > tolerance)
	{
		propagateStep();
		steps += 1;
		conv_prev = conv_curr;
		conv_curr = wvfxn->getNorm();
		wvfxn->normalize();
		conv_ratio = abs(conv_curr - conv_prev);
		if (steps > 1000)
		{
			steps = 0;
			itpdt *= 10.0;
			updateITP();
		}
	}
	copy_n(wvfxn->data(),nx,TestVecs[0]->data());

//Reset variables, get excited states
	if (nv > 1)
	{
		for (int ii = 0; ii < nv; ii++)
		{
			copy_n(TestVecs[ii]->data(),nx,wvfxn->data());
			steps      = 0;
			conv_curr  = 1.0;
			conv_prev  = 1.0;
			conv_ratio = conv_curr/conv_prev;
			while (conv_ratio > tolerance)
			{
				for (int jj = 0; jj < ii; jj++)
				{
					*wvfxn -= (*TestVecs[jj] * (*TestVecs[jj]|*TestVecs[ii]).integrate_rect());
				}
			}
		}
	}

}

void SplitOp1D::initializeGuess(vector<shared_ptr<wvfxn1D>> TestVecs)
{
	int nv = TestVecs.size();
	double span  = 0.5*(xgrid->element(nx-1) - xgrid->element(0));
	double xcen  = 0.5*(xgrid->element(nx-1) + xgrid->element(0));
	double dspan = span/nv;
	double x;
	shared_ptr<Array1D<cplx>> vec;
	for (int ii = 0; ii<nv; ii++)
	{
		vec = TestVecs[ii];
		for (int jj = 0; jj < nx; jj++)
		{
			x = (xgrid->element(jj) - xcen)/(4.0*dspan);
			vec->element(jj) = exp(-1.0*pow(x-0.5,2));
		}
	}
}
