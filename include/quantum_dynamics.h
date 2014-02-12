//	Functions and subroutines for quantum dynamics calculations

typedef std::complex<double> cplx;
using namespace std;

//	Steps a 2D wavefunction forward via a single time step
template <typename T>
void SplitOp2D_Step (Array_2D <T> &A, Array_2D <T> &KinetOp, Array_2D <T> &PotenOp,fftw_plan forplan,fftw_plan backplan)
{
	int ii;
//	Fourier Transform to momentum space, Apply exp(-iTdt/2), invert, normalize
	fftw_execute(forplan);
	for (ii = 0; ii < (A.Nx * A.Ny); ii++) //Arrays are ordered the same so element by element multiplication works
	{
		A.grid[ii] *= KinetOp.grid[ii];
	}
	fftw_execute(backplan);
	fftw_normalize_2D(A);

//	Apply exp(-iVdt)
	for (ii = 0; ii < (A.Nx * A.Ny); ii++)
	{
		A.grid[ii] *= PotenOp.grid[ii];
	}

//	Apply exp(-iTdt/2)
	fftw_execute(forplan);
	for (ii = 0; ii < (A.Nx * A.Ny); ii++)
	{
		A.grid[ii] *= KinetOp.grid[ii];
	}
	fftw_execute(backplan);
	fftw_normalize_2D(A);
	return;
};

//Steps a 1D wavefunction forward via a single time step
template <typename T>
void SplitOp1D_Step (Array_1D <T> &A, Array_1D <T> &KinetOp, Array_1D <T> &PotenOp,fftw_plan forplan,fftw_plan backplan)
{
	int ii;

//	Fourier Transform to momentum space, Apply exp(-iTdt/2), invert, normalize
	fftw_execute(forplan);
	for (ii = 0; ii < A.Nx; ii++) //Arrays are ordered the same so element by element multiplication works
	{
		A.grid[ii] *= KinetOp.grid[ii];
	}
	fftw_execute(backplan);
	fftw_normalize_1D(A);

//	Apply exp(-iVdt)
	for (ii = 0; ii < A.Nx; ii++)
	{
		A.grid[ii] *= PotenOp.grid[ii];
	}

//	Apply exp(-iTdt/2)
	fftw_execute(forplan);
	for (ii = 0; ii < A.Nx; ii++)
	{
		A.grid[ii] *= KinetOp.grid[ii];
	}
	fftw_execute(backplan);
	fftw_normalize_1D(A);
	return;
};

//	Prints the 2D wavefunction to a file, skips some points to minimize file size
template <typename T, typename U>
void print_wvfxn(Array_2D <T> &A, Array_1D <U> &X, Array_1D <U> &Y, int tindex)
{
	ofstream outdata;
	int ii,jj;
    double re,im;

    char filename[sizeof("output_data/wvfxn_t000000.txt")];
    sprintf(filename,"output_data/wvfxn_t%.6d.txt",tindex); //Appended filename
    cout << "filename =  " << filename << endl;
    outdata.open(filename);

    for (ii = 0; ii < A.Nx; ii += 5)
    {
		for (jj = 0; jj < A.Ny; jj += 4)
		{
			re = real(A.get_elem(ii,jj));
			im = imag(A.get_elem(ii,jj));
			outdata << scientific << X.grid[ii] << "\t" << Y.grid[jj] << "\t" << re << "\t" << im << endl;
		}
		outdata << endl;
		}
	outdata.close();
	return;
};

//	Prints a 1D wavefunction slice to a file
template <typename T, typename U>
void print_wvfxn_slice(Array_2D <T> &A, Array_1D <U> &X, int ind, int tindex)
{
	ofstream outdata;
	int ii,jj;
    double re,im;

    char filename[sizeof("output_data/wvfxnslice_t000000.txt")];
    sprintf(filename,"output_data/wvfxnslice_t%.6d.txt",tindex); //Appended filename
    cout << "filename =  " << filename << endl;
    outdata.open(filename);

    for (ii = 0; ii < A.Nx; ii += 5)
    {
		re = real(A.get_elem(ii,ind));
		im = imag(A.get_elem(ii,ind));
		outdata << scientific << X.grid[ii] << "\t" << re << "\t" << im << endl;
	}
	outdata.close();
	return;
};

//	Returns the expectation value in the molecular direction
template <typename T, typename U>
double Z_ExpVal_2D(Array_2D <T> &A, Array_1D <U> &Y, Limits_2D &L, Array_2D <T> &S)
{
	int ii,jj;
	#pragma omp parallel for default(shared) private(ii,jj)
	for (ii = 0; ii < A.Nx; ii++)
	{
		for (jj = 0; jj < A.Ny; jj++)
		{
			S.set_elem(ii,jj,Y.grid[jj]*abs(A.get_elem(ii,jj))*abs(A.get_elem(ii,jj)));
		}
	}
	T zexp = integrate_2D(S,L);
	return real((cplx)zexp);
};

//	Returns the expectation value in the molecular direction
template <typename T, typename U>
double Z_ExpVal_1D(Array_1D <T> &A, Array_1D <U> &Y)
{
	int ii;
	Array_1D <T> wvfxn(A.Nx);
	wvfxn.xstep = A.xstep;
	for (ii = 0; ii < A.Nx; ii++)
	{
		wvfxn.grid[ii] = Y.grid[ii]*abs(A.grid[ii])*abs(A.grid[ii]);
	}
	T zexp = integrate_1D(wvfxn);
	return real((cplx)zexp);
};

//	Calculates a wavefunction correlation function
template <typename T>
cplx correlation(Array_2D <T> &A, Array_2D <T> &B, Array_2D <T> &S)
{
	int ii;
	for (ii = 0; ii < A.Nx*A.Ny; ii++)
    {
    	S.grid[ii] = A.grid[ii] * conj(B.grid[ii]);
    }
    cplx correlation = integrate_2D(S);
	return correlation;
};

//	Overlaps the molecular eigenstates with a cut of the 2D wavefunction in time
template <typename T, typename U>
void get_Dcoefficients(Array_2D <T> &A, EigenArray_1D <T> &Eigen, Array_1D <U> &Storage, int zz)
{
	T * fxn_ptr;
	int ii,jj,kk;
	T dn;
	Array_1D <T> overlap(A.Ny);
	overlap.xstep = A.ystep;
	for (ii = 0; ii < Eigen.n; ii++)
	{
		fxn_ptr = Eigen.get_array_addr(ii);
		for (jj = 0; jj < A.Ny; jj++)
		{
			overlap.grid[jj] = conj((cplx)fxn_ptr[jj]);
		}
		for (kk = 0; kk < A.Ny; kk++)
		{
			overlap.grid[kk] *= A.get_elem(zz,kk);
		}
		dn = integrate_1D(overlap);
		Storage.grid[ii] = dn;
	}
	return;
};


//	Calculates where the energy in the system is currently residing, in eV
template <typename T, typename U>
void get_energy_components(Array_2D <T> &A, Array_1D <U> &X, Array_1D <U> &Y, Array_1D <U> &OutArray, Array_1D <U> &JuncOutArray, Limits_2D &CellLims, Limits_2D &L, Phys_Parameters &P, Array_2D <T> &S)
{
	int ii, jj;
	cplx energy = 0.0;
	T num;
	T poten;

//	Electronic Kinetic Energy
	#pragma omp parallel for default(shared) private(ii,jj,num)
	for (ii = 0; ii < A.Nx; ii++)
	{
		for (jj = 0; jj < A.Ny; jj++)
		{
			num = partial_xx(A,ii,jj);
			S.set_elem(ii,jj,num*conj(A.get_elem(ii,jj)));
		}
	}
	energy = integrate_2D(S,CellLims);
	energy *= (-1.0*P.hbar*P.hbar/(2.0*P.m));
	OutArray.grid[0] = real(energy)*27.21;

	energy = integrate_2D(S,L);
	energy *= (-1.0*P.hbar*P.hbar/(2.0*P.m));
	JuncOutArray.grid[0] = real(energy)*27.21;

//	Molecular Kinetic Energy
	#pragma omp parallel for default(shared) private(ii,jj,num)
	for (ii = 0; ii < A.Nx; ii++)
	{
		for (jj = 0; jj < A.Ny; jj++)
		{
			num = partial_yy(A,ii,jj);
			S.set_elem(ii,jj,num*conj(A.get_elem(ii,jj)));
		}
	}
	energy = integrate_2D(S,CellLims);
	energy *= (-1.0*P.hbar*P.hbar/(2.0*P.M));
	OutArray.grid[1] = real(energy)*27.21;

	energy = integrate_2D(S,L);
	energy *= (-1.0*P.hbar*P.hbar/(2.0*P.M));
	JuncOutArray.grid[1] = real(energy)*27.21;

//	Electronic Potential Energy
	#pragma omp parallel for default(shared) private(ii,jj,num,poten)
	for (ii = 0; ii < A.Nx; ii++)
	{
		poten = V_e(X.grid[ii],P);
		for (jj = 0; jj < A.Ny; jj++)
		{
			num = A.get_elem(ii,jj);
			S.set_elem(ii,jj,num*conj(num)*poten);
		}
	}
	energy = integrate_2D(S,CellLims);
	OutArray.grid[2] = real(energy)*27.21;

	energy = integrate_2D(S,L);
	JuncOutArray.grid[2] = real(energy)*27.21;


//	Molecular Potential Energy
	#pragma omp parallel for default(shared) private(ii,jj,num,poten)
	for (jj = 0; jj < A.Ny; jj++)
	{
		poten = V_M(Y.grid[jj],P);
		for (ii = 0; ii < A.Nx; ii++)
		{
			num = A.get_elem(ii,jj);
			S.set_elem(ii,jj,num*conj(num)*poten);
		}
	}
	energy = integrate_2D(S,CellLims);
	OutArray.grid[3] = real(energy)*27.21;

	energy = integrate_2D(S,L);
	JuncOutArray.grid[3] = real(energy)*27.21;

//	Interaction Energy
	#pragma omp parallel for default(shared) private(ii,jj,num,poten)
	for (ii = 0; ii < A.Nx; ii++)
	{
		for (jj = 0; jj < A.Ny; jj++)
		{
			poten = V_eM(X.grid[ii],Y.grid[jj],P);
			num = A.get_elem(ii,jj);
			S.set_elem(ii,jj,num*conj(num)*poten);
		}
	}
	energy = integrate_2D(S,CellLims);
	OutArray.grid[4] = real(energy)*27.21;

	energy = integrate_2D(S,L);
	JuncOutArray.grid[4] = real(energy)*27.21;

	return;
};

//	Calculates the wavefunction probability density in the junction region
template <typename T>
double get_junction_wavefunction(Array_2D <T> &A, Limits_2D &L, Array_2D <T> &S)
{
	int ii,jj;
	#pragma omp parallel for default(shared) private(ii,jj)
	for (ii = 0; ii < A.Nx; ii++)
	{
		for (jj = 0; jj < A.Ny; jj++)
		{
			S.set_elem(ii,jj,conj(A.get_elem(ii,jj))*A.get_elem(ii,jj));
		}
	}

	cplx prob_density 	= integrate_2D(S,L);
	double re 			= real(prob_density);
//	cout << "The Probability density in the junction region is " << re << endl;
	return re;
};

//	Calculate the C_nm coefficients within the junction region
//	Offset term takes into account that the first value in ElecStates and Psi correspond to different x values
template <typename T>
cplx get_junction_Cnm(Array_2D <T> &Psi, EigenArray_1D <T> &ElecStates, EigenArray_1D <T> &MolStates, Limits_2D &L, int m, int n, int offset, Array_2D <T> &Scrap)
{
	int ii,jj;
	T *Elec 	= ElecStates.get_array_addr(n);  
	T *Mol 		= MolStates.get_array_addr(m);

	#pragma omp parallel for default(shared) private(ii,jj)
	for (ii = 0; ii < ElecStates.Nx; ii++)
	{
		for (jj = 0; jj < MolStates.Nx; jj++)
		{
			Scrap.set_elem(ii+offset,jj, Psi.get_elem(ii+offset,jj)*Elec[ii]*Mol[jj]); 
		}
	}

	cplx Cnm 	= integrate_2D(Scrap,L); 
	return Cnm;
}
