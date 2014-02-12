//	Set of numeric algorithms for analysis of members of the Array_nD class, other mathematical functions

typedef std::complex<double> cplx;
using namespace std;

//	Sixth order accurate numeric first derivative, 1D array
template <typename T>
T deriv_1D_1st(Array_1D <T> &A, int N)
{
	if (N < 3) 
	{
		T deriv = (1/(12*A.xstep)) * (-25*A.grid[N] + 48*A.grid[N+1]-36*A.grid[N+2]+16*A.grid[N+3]-3*A.grid[N+4]);
		return deriv;
	} 
	else if (N > (A.Nx-4))
	{
		T deriv = (1/(12*A.xstep)) * (25*A.grid[N]-48*A.grid[N-1]+36*A.grid[N-2]-16*A.grid[N-3]+3*A.grid[N-4]);
		return deriv;
	} 
	else 
	{
		T deriv = (1.0/(60*A.xstep)) * (-1.0*A.grid[N-3] + 9.0*A.grid[N-2] - 45.0*A.grid[N-1] + 45.0*A.grid[N+1] - 9.0*A.grid[N+2] + 1.0*A.grid[N+3]);
		return deriv;
	}
};

//	Sixth order accurate numeric second derivative, 1D array
template <typename T>
T deriv_1D_2nd(Array_1D <T> &A, int N)
{
	if (N < 3) 
	{
		T deriv = (45*A.grid[N] - 154*A.grid[N+1] + 214*A.grid[N+2] - 156*A.grid[N+3] + 61*A.grid[N+4] - 10*A.grid[N+5]);
		deriv /= (12*A.xstep*A.xstep);
		return deriv;
	}
	else if (N > (A.Nx-4)) 
	{
		T deriv = (45*A.grid[N] - 154*A.grid[N-1] + 214*A.grid[N-2] - 156*A.grid[N-3] + 61*A.grid[N-4] - 10*A.grid[N-5]);
		deriv /= (12*A.xstep*A.xstep);
		return deriv;
	}
	else 
	{
		T deriv = 2*A.grid[N-3] - 27*A.grid[N-2] + 270*A.grid[N-1] - 490*A.grid[N] + 270*A.grid[N+1] - 27*A.grid[N+2] + 2*A.grid[N+3];
		deriv /= (180*A.xstep*A.xstep);
		return deriv;
	}
};

//	Numeric (trapezoidal) integration, 1D array, provide limits
template <typename T>
T integrate_1D(Array_1D <T> &A, int NLower, int NUpper)
{
	if (NUpper >= A.Nx) 
	{
		NUpper = A.Nx-1; //correct if upper limit is too large
	}
	if (NLower < 0 && NLower >= NUpper) 
	{
		NLower = 0; //correct if lower limit is too large or negative
	}
	T sum = 0;
	int ii;
	for (ii = NLower+1; ii < NUpper; ii++) 
	{
		sum += A.grid[ii];
	}
	sum += (A.grid[NLower] + A.grid[NUpper])/2;
	sum *= A.xstep;
	return sum;
};

//	Numeric (trapezoidal) integration, 1D array, full array integration
template <typename T>
T integrate_1D(Array_1D <T> &A)
{
	T sum = 0;
	int ii;
	for (ii = 1; ii < A.Nx-1; ii++) 
	{
		sum += A.grid[ii];
	}
	sum += (A.grid[0] + A.grid[A.Nx-1])/2;
	sum *= A.xstep;
	return sum;
};

//	Sixth order accurate numeric first derivative, 2D array, x-direction
template <typename T>
T partial_x(Array_2D <T> &A, int Nx, int Ny)
{
	if (Nx < 3) 
	{
		T deriv = (1/(12*A.xstep)) * (-25*A.get_elem(Nx,Ny) + 48*A.get_elem(Nx+1,Ny) - 36*A.get_elem(Nx+2,Ny) + 16*A.get_elem(Nx+3,Ny) - 3*A.get_elem(Nx+4,Ny));
		return deriv;
	}
	else if (Nx > (A.Nx-4))
	{
		T deriv = (1/(12*A.xstep)) * (25*A.get_elem(Nx,Ny) - 48*A.get_elem(Nx-1,Ny) + 36*A.get_elem(Nx-2,Ny) - 16*A.get_elem(Nx-3,Ny) + 3*A.get_elem(Nx-4,Ny));
		return deriv;
	}
	else 
	{
		T deriv = (1.0/(60*A.xstep)) * (-1.0*A.get_elem(Nx-3,Ny) + 9.0*A.get_elem(Nx-2,Ny) - 45.0*A.get_elem(Nx-1,Ny) + 45.0*A.get_elem(Nx+1,Ny) - 9.0*A.get_elem(Nx+2,Ny) + 1.0*A.get_elem(Nx+3,Ny));
		return deriv;
	}
};

//	Sixth order accurate numeric first derivative, 2D array, y-direction
template <typename T>
T partial_y(Array_2D <T> &A, int Nx, int Ny)
{
	if (Ny < 3) 
	{
		T deriv = (1/(12.0*A.ystep)) * (-25*A.get_elem(Nx,Ny) + 48*A.get_elem(Nx,Ny+1) - 36*A.get_elem(Nx,Ny+2) + 16*A.get_elem(Nx,Ny+3) - 3*A.get_elem(Nx,Ny+4));
		return deriv;
	}
	else if (Ny > (A.Ny-4))
	{
		T deriv = (1/(12.0*A.ystep)) * (25*A.get_elem(Nx,Ny) - 48*A.get_elem(Nx,Ny-1) + 36*A.get_elem(Nx,Ny-2) - 16*A.get_elem(Nx,Ny-3) + 3*A.get_elem(Nx,Ny-4));
		return deriv;
	}
	else 
	{
		T deriv = (1.0/(60.0*A.ystep)) * (-1.0*A.get_elem(Nx,Ny-3) + 9.0*A.get_elem(Nx,Ny-2) - 45.0*A.get_elem(Nx,Ny-1) + 45.0*A.get_elem(Nx,Ny+1) - 9.0*A.get_elem(Nx,Ny+2) + 1.0*A.get_elem(Nx,Ny+3));
		return deriv;
	}
};

//	Fourth order accurate numeric second derivative, 2D array, x-direction
template <typename T>
T partial_xx(Array_2D <T> &A, int Nx, int Ny)
{
	if (Nx < 3) 
	{
		T deriv = (1/(A.xstep*A.xstep)) * (\
			(15.0/4.0)*A.get_elem(Nx,Ny) \
			- (77.0/6.0)*A.get_elem(Nx+1,Ny) \
			+ (107.0/6.0)*A.get_elem(Nx+2,Ny) \
			- (13.0)*A.get_elem(Nx+3,Ny) \
			+ (61.0/12.0)*A.get_elem(Nx+4,Ny)\
			- (5.0/6.0)*A.get_elem(Nx+5,Ny));
		return deriv;
	}
	else if (Nx > (A.Nx-4))
	{
		T deriv = (1/(A.xstep*A.xstep)) * (\
			(15.0/4.0)*A.get_elem(Nx,Ny) \
			- (77.0/6.0)*A.get_elem(Nx-1,Ny) \
			+ (107.0/6.0)*A.get_elem(Nx-2,Ny) \
			- (13.0)*A.get_elem(Nx-3,Ny) \
			+ (61.0/12.0)*A.get_elem(Nx-4,Ny) \
			- (5.0/6.0)*A.get_elem(Nx-5,Ny));
		return deriv;
	}
	else 
	{
		T deriv = (1.0/(A.xstep*A.xstep)) * (\
			(1.0/90.0)*A.get_elem(Nx-3,Ny) \
			- (3.0/20.0)*A.get_elem(Nx-2,Ny) \
				+ (3.0/2.0)*A.get_elem(Nx-1,Ny) \
			- (49.0/18.0)*A.get_elem(Nx,Ny) \
			+ (3.0/2.0)*A.get_elem(Nx+1,Ny) \
			- (3.0/20.0)*A.get_elem(Nx+2,Ny) \
			+ (1.0/90.0)*A.get_elem(Nx+3,Ny));
		return deriv;
	}
};

//	Second order accurate numeric second derivative, 2D array, x-direction
template <typename T>
T partial_xx_2(Array_2D <T> &A, int Nx, int Ny)
{
	if (Nx < 3) 
	{
		T deriv = (1/(A.xstep*A.xstep)) * (\
			2.0*A.get_elem(Nx,Ny) \
			- 5.0*A.get_elem(Nx+1,Ny) \
			+ 4.0*A.get_elem(Nx+2,Ny) \
			- 1.0*A.get_elem(Nx+3,Ny));
		return deriv;
	} 
	else if (Nx > (A.Nx-4))
	{
		T deriv = (1/(A.xstep*A.xstep)) * (\
			2.0*A.get_elem(Nx,Ny) \
			- 5.0*A.get_elem(Nx-1,Ny) \
			+ 4.0*A.get_elem(Nx-2,Ny) \
			- 1.0*A.get_elem(Nx-3,Ny));
		return deriv;
	}
	else 
	{
		T deriv = (1.0/(A.xstep*A.xstep)) * (\
			- (1.0/12.0)*A.get_elem(Nx-2,Ny) \
			+ (4.0/3.0)*A.get_elem(Nx-1,Ny) \
			- (5.0/2.0)*A.get_elem(Nx,Ny) \
			+ (4.0/3.0)*A.get_elem(Nx+1,Ny) \
			- (1.0/12.0)*A.get_elem(Nx+2,Ny));
		return deriv;
	}
};

//	Fourth order accurate numeric second derivative, 2D array, y-direction
template <typename T>
T partial_yy(Array_2D <T> &A, int Nx, int Ny)
{
	if (Ny < 3) 
	{
		T deriv = (1/(A.ystep*A.ystep)) * (\
			(15.0/4.0)*A.get_elem(Nx,Ny) \
			- (77.0/6.0)*A.get_elem(Nx,Ny+1) \
			+ (107.0/6.0)*A.get_elem(Nx,Ny+2) \
			- (13.0)*A.get_elem(Nx,Ny+3) \
			+ (61.0/12.0)*A.get_elem(Nx,Ny+4)\
			- (5.0/6.0)*A.get_elem(Nx,Ny+5));
		return deriv;
	} 
	else if (Ny > (A.Ny-4))
	{
		T deriv = (1/(A.ystep*A.ystep)) * (\
			(15.0/4.0)*A.get_elem(Nx,Ny) \
			- (77.0/6.0)*A.get_elem(Nx,Ny-1) \
			+ (107.0/6.0)*A.get_elem(Nx,Ny-2) \
			- (13.0)*A.get_elem(Nx,Ny-3) \
			+ (61.0/12.0)*A.get_elem(Nx,Ny-4)\
			- (5.0/6.0)*A.get_elem(Nx,Ny-5));
		return deriv;
	} 
	else 
	{
		T deriv = (1.0/(A.ystep*A.ystep)) * (\
			(1.0/90.0)*A.get_elem(Nx,Ny-3) \
			- (3.0/20.0)*A.get_elem(Nx,Ny-2) \
			+ (3.0/2.0)*A.get_elem(Nx,Ny-1) \
			- (49.0/18.0)*A.get_elem(Nx,Ny) \
			+ (3.0/2.0)*A.get_elem(Nx,Ny+1) \
			- (3.0/20.0)*A.get_elem(Nx,Ny+2) \
			+ (1.0/90.0)*A.get_elem(Nx,Ny+3));
		return deriv;
	}
};

//	Numeric (trapezoidal) integration, 2D array, full array integration
//	1 2 2 2 1
//	2 4 4 4 2
//  2 4 4 4 2
//  1 2 2 2 1
template <typename T>
T integrate_2D(Array_2D <T> &A)
{
	int NX 	= A.Nx-1;
	int NY 	= A.Ny-1;
	T sum1 	= 0;
	T sum2 	= 0;
	T sum3 	= 0;
	int ii,jj;

//	Corner terms
	sum1 = (A.get_elem(0,0) + A.get_elem(0,NY) + A.get_elem(NX,0) + A.get_elem(NX,NY));

//	Border terms
	for (ii = 1; ii < NX; ii++) 
	{
		sum2 += A.get_elem(ii,0);
	}
	for (ii = 1; ii < NX; ii++) 
	{
		sum2 += A.get_elem(ii,NY);
	}
	for (jj = 1; jj < NY; jj++) 
	{
		sum2 += A.get_elem(0,jj);
	}
	for (jj = 1; jj < NY; jj++) 
	{
		sum2 += A.get_elem(NX,jj);
	}
	sum2 *= 2;

//	Center terms
	for (ii = 1; ii<NX; ii++)
	{
		for (jj = 1; jj<NY; jj++) 
		{
			sum3 += A.get_elem(ii,jj);
		}
	}
	sum3 *= 4;

//	Finalize
	T sum = A.xstep * A.ystep * 0.25 * (sum1 + sum2 + sum3);
	return sum;
};

//	Numeric (trapezoidal) integration, 2D array, provide limits
template <typename T>
T integrate_2D(Array_2D <T> &A, Limits_2D &L)
{
	T sum1 = 0;
	T sum2 = 0;
	T sum3 = 0;
	int ii,jj;

//	Corner terms
	sum1 = (A.get_elem(L.XL,L.YL) + A.get_elem(L.XL,L.YU) + A.get_elem(L.XU,L.YL) + A.get_elem(L.XU,L.YU));

//	Border terms
	for (ii = L.XL+1; ii < L.XU; ii++) 
	{
		sum2 += A.get_elem(ii,L.YL);
	}
	for (ii = L.XL+1; ii < L.XU; ii++) 
	{
		sum2 += A.get_elem(ii,L.YU);
	}
	for (jj = L.YL+1; jj < L.YU; jj++) 
	{
		sum2 += A.get_elem(L.XL,jj);
	}
	for (jj = L.YL+1; jj < L.YU; jj++) 
	{
		sum2 += A.get_elem(L.XU,jj);
	}
	sum2 *= 2;

//	Center terms
	for (ii = L.XL+1; ii<L.XU; ii++)
	{
		for (jj = L.YL+1; jj<L.YU; jj++) 
		{
			sum3 += A.get_elem(ii,jj);
		}
	}
	sum3 *= 4;

//	Finalize
	T sum = A.xstep * A.ystep * 0.25 * (sum1 + sum2 + sum3);
	return sum;
};

//	Probability current for a 2D array, x direction
template <typename T>
double probcurrent_x_2D(Array_2D <T> &A, int xx, int yy, double mass, double hbar)
{
	cplx I 			= 1.0j;
	cplx psi_star 	= conj((cplx)A.get_elem(xx,yy));
	cplx psi_deriv 	= -1.0 * I * (cplx)partial_x(A,xx,yy);
	double probcur 	= real(psi_star * psi_deriv)*hbar/mass;
	return probcur;
};

//	Probability current for a 2D array, y direction
template <typename T>
double probcurrent_y_2D(Array_2D <T> &A, int xx, int yy, double mass, double hbar)
{
	cplx I 			= 1.0j;
	cplx psi_star 	= conj((cplx)A.get_elem(xx,yy));
	cplx psi_deriv 	= -1.0 * I * (cplx)partial_y(A,xx,yy);
	double probcur 	= real(psi_star * psi_deriv)*hbar/mass;
	return probcur;
};

//	Probability current for a 1D array, used in flux calculations
template <typename T>
double probcurrent_1D(Array_1D <T> &A, int xx, double mass, double hbar)
{
	cplx I 			= 1.0j;
	cplx psi_star 	= conj((cplx)A.grid[xx]);
	cplx psi_deriv 	= -1.0 * I * (cplx)deriv_1D_1st(A,xx);
	double probcur 	= real(psi_star * psi_deriv)*hbar/mass;
	return probcur;
};

//	Probability current for a 2D array, used in flux calculations
template <typename T>
double probcurrent_2D(Array_2D <T> &A, int xx, double mass, double hbar)
{
	int ii;
	double sum 	= 0.0;
	cplx I 		= 1.0j;
	cplx psi_star, psi_deriv;
	double probcur;
	Array_1D <double> probcurrent(A.Ny);
	probcurrent.xstep = A.ystep;
	for (ii = 0; ii < A.Ny; ii++)
	{
		psi_star 				= conj((cplx)A.get_elem(xx,ii));
		psi_deriv 				= -1.0 * I * (cplx)partial_x(A,xx,ii);
		probcur 				= real(psi_star * psi_deriv)*hbar/mass;
		probcurrent.grid[ii] 	= probcur;
	}
	probcur = integrate_1D(probcurrent);
	return probcur;
};

//	Returns the nth hermite polynomial evaluated at x
double hermite(int,double);

//	Factorial function
double factorial(int);

//	Returns a harmonic oscillator wavefunction with the given parameters at a given point in space
double ho_wvfxn(int, double, double, double, double, double); //Harmonic oscillator wavefunctions

//	Kronecker Delta function
int kron_delta(int, int);

//	Calculate the value of \psi^* \psi integrate over all space, 2D
template <typename T>
double get_wvfxn_norm_2D(Array_2D <T> &A)
{
	int ii;
	double Re, Im;
	Array_2D <double> wvfxn_norm(A.Nx,A.Ny);
	wvfxn_norm.xstep = A.xstep;
	wvfxn_norm.ystep = A.ystep;
	for (ii = 0; ii < A.Nx*A.Ny; ii++)
	{
		Re 					= real((cplx)A.grid[ii]);
		Im 					= imag((cplx)A.grid[ii]);
		wvfxn_norm.grid[ii] = Re*Re + Im*Im;
	}
	double norm = sqrt(integrate_2D(wvfxn_norm));
	return norm;
};

//	Same as previous function, but with additional argument for pre-allocated array
template <typename T>
double get_wvfxn_norm_2D(Array_2D <T> &A, Array_2D <T> &S)
{
	int ii;
	double Re, Im;
	for (ii = 0; ii < A.Nx*A.Ny; ii++){
		Re 			= real((cplx)A.grid[ii]);
		Im 			= imag((cplx)A.grid[ii]);
		S.grid[ii] 	= Re*Re + Im*Im;
	}
	cplx num 		= integrate_2D(S);
	double norm 	= sqrt(num.real());
	return norm;
};

//	Same as previous function, but with additional argument to provide limits of integration
template <typename T>
double get_wvfxn_norm_2D(Array_2D <T> &A, Limits_2D &L, Array_2D <T> &S)
{
	int ii;
	double Re, Im;
	#pragma omp parallel for default(shared) private(Re,Im,ii)
	for (ii = 0; ii < A.Nx*A.Ny; ii++){
		Re 			= real((cplx)A.grid[ii]);
		Im 			= imag((cplx)A.grid[ii]);
		S.grid[ii] 	= Re*Re + Im*Im;
	}
	cplx num 		= integrate_2D(S,L);
	double norm 	= sqrt(num.real());
	return norm;
};

//	Calculate the value of \psi^* \psi integrate over all space, 1D
template <typename T>
double get_wvfxn_norm_1D(Array_1D <T> &A)
{
	int ii;
	double Re, Im;
	Array_1D <double> wvfxn_norm(A.Nx);
	wvfxn_norm.xstep = A.xstep;
	for (ii = 0; ii < A.Nx; ii++)
	{
		Re 					= real((cplx)A.grid[ii]);
		Im 					= imag((cplx)A.grid[ii]);
		wvfxn_norm.grid[ii] = Re*Re + Im*Im;
	}
	double norm = integrate_1D(wvfxn_norm);
	norm 		= sqrt(norm);
	return norm;
};

//	Normalize a 1D wavefunction to 1.00
template <typename T>
void normalize_wxfxn_1D(Array_1D <T> &A)
{
	int ii;
	double norm = get_wvfxn_norm_1D(A);
	for (ii = 0; ii < A.Nx; ii++)
	{
		A.grid[ii] *= 1.0/norm;
	}
	return;
};

//	Normalize a 2D wavefunction to 1.00
template <typename T>
void normalize_wxfxn_2D(Array_2D <T> &A)
{
	int ii;
	double norm = get_wvfxn_norm_2D(A);
	for (ii = 0; ii < A.Nx*A.Ny; ii++) 
	{
		A.grid[ii] /= norm;
	}
};

//	Same as previous function, but with additional argument for pre-allocated array
template <typename T>
void normalize_wxfxn_2D(Array_2D <T> &A, Array_2D <T> &S)
{
	int ii;
	double norm = get_wvfxn_norm_2D(A,S);
	for (ii = 0; ii < A.Nx*A.Ny; ii++) 
	{
		A.grid[ii] /= norm;
	}
};
