typedef std::complex<double> cplx;
using namespace std;
using namespace boost;

//	Single imaginary time propagation step on 1D array;
template <typename T, typename U>
void ITP_single_step(double dt, U *ptr, Array_1D <T> &PotenOp, Array_1D <T> &KinetOp)
{
//	Put on a single processor for a 1D calculation
	int ii;

//	Set up fourier transforms
	fftw_plan_with_nthreads(1);
	fftw_plan forplan, backplan;
	forplan 	= fftw_plan_dft_1d(PotenOp.Nx,(fftw_complex*)ptr,(fftw_complex*)ptr,FFTW_FORWARD,FFTW_ESTIMATE);
	backplan 	= fftw_plan_dft_1d(PotenOp.Nx,(fftw_complex*)ptr,(fftw_complex*)ptr,FFTW_BACKWARD,FFTW_ESTIMATE);

//	Step in imaginary time

//	Momentum space
	fftw_execute(forplan);
	for (ii = 0; ii < PotenOp.Nx; ii++)
	{
		ptr[ii] *= KinetOp.grid[ii];
	}
	fftw_execute(backplan);
	fftw_normalize_1D(ptr,PotenOp.Nx);

//	Position Space
	for (ii = 0; ii < PotenOp.Nx; ii++)
	{
		ptr[ii] *= PotenOp.grid[ii];
	}

//	Momentum space
	fftw_execute(forplan);
	for (ii = 0; ii < PotenOp.Nx; ii++)
	{
		ptr[ii] *= KinetOp.grid[ii];
	}
	fftw_execute(backplan);
	fftw_normalize_1D(ptr,PotenOp.Nx);

//	Clean Up
	fftw_destroy_plan(forplan);
	fftw_destroy_plan(backplan);

	return;
};

//	Fills a 1D EigenArray with initial function, set to be Gaussian curves
template <typename U, typename T>
void initialize_fxn_array(EigenArray_1D <U> &A, Array_1D <T> &xgrid)
{
	int ii,jj,kk;
	T span 	= (xgrid.grid[xgrid.Nx-1] - xgrid.grid[0]) / 2.0;
	T xcen 	= (xgrid.grid[xgrid.Nx-1] + xgrid.grid[0]) / 2.0;
	T dspan = span / A.n;
	T x;
	U * array_ptr;
	for (ii = 0; ii<A.n; ii++){
		array_ptr = A.get_array_addr(ii);
		for (jj = 0; jj < A.Nx; jj++)
		{
			x 				= (xgrid.grid[jj] - xcen)/(0.4*dspan);
			array_ptr[jj] 	= exp(-1.0*(x-0.5)*(x-0.5));
		}
	}
	A.normalize();

//	Print initialized array to file to check
	ofstream outfile;
	outfile.open("output_data/InitialWavefunctionArray.txt");
	for (jj = 0; jj < A.Nx; jj++)
	{
		outfile << xgrid.grid[jj] << " ";
		for (kk = 0; kk < A.n; kk++)
		{
			outfile << real((cplx)(A.get_array_addr(kk))[jj]) << " ";
		}
		outfile << endl;
	}
	outfile.close();
	return;
  };

//	Projects the n1th array of the EigenArray onto the n2th array
template <typename T>
cplx projection (EigenArray_1D <T> &A, int n1, int n2)
{
	Array_1D <cplx> wvfxn_prod(A.Nx);
	wvfxn_prod.xstep = real((cplx)A.xstep);
	int ii;
	T *u = A.get_array_addr(n1);
	T *v = A.get_array_addr(n2);
	for (ii = 0; ii < A.Nx; ii++)
	{
		wvfxn_prod.grid[ii] = u[ii] * conj((cplx)(v[ii]));
	}
	cplx proj = integrate_1D(wvfxn_prod);
	return proj;
};

//	Load EigenArray from a file, skips first column, loads other available columns
//	Possible to maintain functionality if there are too many lines or columns in file
template <typename T>
void load_eigen_array(EigenArray_1D <T> &A)
{
	int ii, jj;
	jj = 0;
	T * array;
	string line;
	vector <string> strs;
	ifstream infile("output_data/InitialWavefunctionArray.txt");
	if (infile.is_open())
	{
		while (getline(infile,line) && jj < A.Nx)
		{
			split(strs,line,is_any_of("\t \n"));
			for (ii = 1; ii < A.n; ii++) //Starts at 1 since the first column is skipped
			{
				array 		= A.get_array_addr(ii);
				array[jj] 	= atof(strs[ii].c_str());
			}
			jj++;
		}
	}
	else
	{
		cout << "Unable to open eigenarray file" << endl;
		exit(-1);
	}
	infile.close();
	return;
};

//	Calculates the energies of the eigenfunctions determined via imaginary time propagation, prints to screen
//	The potential energy surfaces must be real of type double
//	Uses derivatives rather than the Fourier transform method
template <typename T>
void calculate_energies (EigenArray_1D <T> &A, Array_1D <double> &X, Array_1D <double> &V, Phys_Parameters &P, double mass)
{
	int ii, jj;
	T *ptr;
	double Tsum, Vsum, energy;
	Array_1D <cplx> wvfxn_kinetic(A.Nx);
	wvfxn_kinetic.xstep = X.xstep;
	Array_1D <cplx> wvfxn_kinetic_deriv(A.Nx);
	wvfxn_kinetic_deriv.xstep = X.xstep;
	Array_1D <cplx> wvfxn_potential(A.Nx);
	wvfxn_potential.xstep = X.xstep;
	double E1, E2;
	for (ii = 0; ii < A.n; ii++)
	{
		ptr = A.get_array_addr(ii);
		for (jj = 0; jj < A.Nx; jj++)
		{
			wvfxn_kinetic.grid[jj] 			= real((cplx)ptr[jj]);
			wvfxn_potential.grid[jj] 		= conj((cplx)ptr[jj])*(cplx)ptr[jj]*V.grid[jj];
		}
		for (jj = 0; jj < A.Nx; jj++)
		{
			wvfxn_kinetic_deriv.grid[jj] 	= deriv_1D_2nd(wvfxn_kinetic,jj);
		}
		for (jj = 0; jj < A.Nx; jj++)
		{
			wvfxn_kinetic_deriv.grid[jj] 	*= wvfxn_kinetic.grid[jj] * (-1.0*P.hbar*P.hbar/(2.0*mass));
		}
		Tsum 	= real(integrate_1D(wvfxn_kinetic_deriv));
		Vsum 	= real(integrate_1D(wvfxn_potential));
		energy 	= Vsum + Tsum;
		cout << "Energy of state " << ii << " is " << energy << " Hartree or " << 27.21*energy << " eV." << endl;

		if (A.n > 0)
		{
			if (ii == 0)
			{
				E1 = energy;
			}
			else if (ii == 1)
			{
				E2 = energy;
			}
		}
	}

	//Finds hbar*omega using the difference between the first two levels, approximating to harmonic oscillator solution
	if (A.n > 0)
	{
		cout << "hw = " << abs(E1-E2) << " Hartree (" << 27.21*abs(E1-E2) << " eV)" << endl;
	}

	return;
};

//	A should be of type cplx or fftw_complex
//	T is double or float, U is cplx or fftw_complex
template <typename U>
void get_wvfxns (EigenArray_1D <U> &A, Array_1D <double> &X, Array_1D <double> &V, Array_1D <double> &K, Phys_Parameters &P, double mass)
{
	U *ptr, *proj_ptr, proj;
	int ii,jj,kk,steps;

//	Constants needed for the adaptive time stepping
	double conv_factor 	= 1e-10; //The integral of the wavefunction after an imaginary time step must constant
	double conv_prev 	= 1.0;
	double conv_curr 	= 1.0;
	double conv_ratio 	= conv_curr/conv_prev;
	double dt;

//	Make initial guess for the time step needed
	if (mass < 10.0) {dt = 1e-2;} else {dt = 5.0;}
	double coeff,time;

//	initialize the function array with arbitrarily placed gaussian curves
	initialize_fxn_array(A,X);
	A.normalize();

//	Set up arrays for transformation
	Array_1D <double> KinetOp(A.Nx);
	Array_1D <double> PotenOp(A.Nx);
	for (ii = 0; ii < A.Nx; ii++)
	{
		KinetOp.grid[ii] = exp(-0.5*dt*K.grid[ii]);
		PotenOp.grid[ii] = exp(-1.0*dt*V.grid[ii]);
	}

//	Find the ground state wavefunction
	cout << "Calculating the ground state..." << endl;
	steps 	= 0;
	ptr 	= A.get_array_addr(0);
	while (conv_ratio > conv_factor)
	{
//	Single time step
		ITP_single_step(dt,ptr,PotenOp,KinetOp);
		steps++;
		conv_prev = conv_curr;
		conv_curr = A.normalize_n(0); //This step normalizes the wavefunction
//		cout << "\r" << conv_curr << flush;
		conv_ratio = abs(conv_curr - conv_prev); //Print to file?

//	Procedure for adjusting the time step
		if (steps > 1000)
		{
			steps = 0;
			dt *= 10.0;
			for (ii = 0; ii < A.Nx; ii++)
			{
				KinetOp.grid[ii] = exp(-0.5*dt*K.grid[ii]);
				PotenOp.grid[ii] = exp(-1.0*dt*V.grid[ii]);
			}
		}

	}
	cout << endl;

//	Reset values
	steps = 0;
	conv_prev = 1.0;
	conv_curr = 1.0;
	conv_ratio = conv_curr/conv_prev;
	if (mass < 10.0) {dt = 1e-2;} else {dt = 5.0;}

//	Find the excited state wavefunctions

	if (A.n > 1)
	{
		for (ii = 1; ii < A.n; ii++)
		{
			steps 		= 0;
			conv_prev 	= 1.0;
			conv_curr 	= 1.0;
			conv_ratio 	= conv_curr/conv_prev;
			if (mass < 10.0) {dt = 1e-2;} else {dt = 5.0;}

			ptr = A.get_array_addr(ii);
			cout << "Calculating excited state #" << ii << endl;
			while (conv_ratio > conv_factor)
			{
//	Orthonormalization Procedure
				for (jj = 0; jj < ii; jj++)
				{
					proj_ptr 	= A.get_array_addr(jj);
					proj 		= projection(A,ii,jj);
					for (kk = 0; kk < A.Nx; kk++)
					{
						ptr[kk] -= (proj * proj_ptr[kk]);
					}
				}
//	Single time step
				ITP_single_step(dt,ptr,PotenOp,KinetOp);
				steps++;
				conv_prev 	= conv_curr;
				conv_curr 	= A.normalize_n(ii); //This step normalizes the wavefunction
//				cout << "\r" << conv_curr << flush;
				conv_ratio 	= abs(conv_curr - conv_prev); //Print to file?
//	Procedure for adjusting the time step
				if (steps > 1000)
				{
					steps = 0;
					dt *= 2.0;
					for (jj = 0; jj < A.Nx; jj++)
					{
						KinetOp.grid[jj] = exp(-0.5*dt*K.grid[jj]);
						PotenOp.grid[jj] = exp(-1.0*dt*V.grid[jj]);
					}
				}

			}

			cout << endl;
		}
	}
	calculate_energies(A,X,V,P,mass);
	return;
};
