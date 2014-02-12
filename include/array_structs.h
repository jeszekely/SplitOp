typedef std::complex<double> cplx;
using namespace std;

//	Class for a 2D array and associated variables needed for split operator calculations
template <typename T> class Array_2D
{
	public:
		int Nx, Ny;

		double xstep, ystep;
		double xinit, yinit;
		double time;
		T * grid;

		//Constructor
		Array_2D (int xx, int yy)
		{
			Nx 		= xx;
			Ny 		= yy;
			xstep 	= ystep = xinit = yinit = 0.0;
			grid 	= new T [xx*yy];
			time 	= 0.0;
		}

		//Other Constructor
		Array_2D (int xx, int yy, double * vals)
		{
			Nx 		= xx;
			Ny 		= yy;
			xstep 	= vals[0];
			xinit 	= vals[1];
			ystep 	= vals[2];
			yinit 	= vals[3];
			grid 	= new T [xx*yy];
			time 	= 0.0;
		}

		//Destructor
		~Array_2D () {delete [] grid;}

		T get_elem(int kk, int ll) {return grid[kk*Ny + ll];} //return a particular element of the array

		void set_elem(int ii, int jj, T num) {grid[ii*Ny + jj] = num;} //return a particular element of the array

		void set_time(double tt) {time = tt;} //manually set the time variable if needed

};

//	Class for a 1D array and associated variables needed for split operator calculations

template <typename T> class Array_1D
{
	public:
		int Nx;

		double xstep, xinit;
		double time;
		T *grid;

		//Constructor
		Array_1D (int xx)
		{
			Nx 		= xx;
			xstep 	= xinit = 0.0;
			grid 	= new T [xx];
			time 	= 0.0;
		}

		//Other Constructor
		Array_1D (int xx, double * vals)
		{
			Nx 		= xx;
			xstep 	= vals[0];
			xinit 	= vals[1];
			grid 	= new T [xx];
			time 	= 0.0;
		}

		//Destructor
		~Array_1D ()
		{
			delete [] grid;
		}

		T get_elem(int ii) {return grid[ii];} //return a particular element of the array

		void set_elem(int ii, T num) {grid[ii] = num;} //return a particular element of the array

		void set_time(double tt) {time = tt;} //manually set the time variable if needed

		void fill_array()
		{
			int ii;
			for (ii = 0; ii<Nx; ii++) 
			{
				grid[ii] = (xinit + ii*xstep);
			}
		}
};

//	Correct the multiplication done by FFTW, 1D array
template <typename T> void fftw_normalize_1D(Array_1D <T> &A)
{
	int ii;
	int N = A.Nx;
	for (ii=0; ii<N; ii++)
	{
		A.grid[ii] /= N;
	}
}

//	Correct the multiplication done by FFTW, 2D array
template <typename T> void fftw_normalize_2D(Array_2D <T> &A)
{
	int ii;
	int N = A.Nx * A.Ny;
	for (ii=0; ii<N; ii++)
	{
		A.grid[ii] /= N;
	}
}

//	Correct the multiplication done by FFTW, 1D array, normal array class
template <typename T> void fftw_normalize_1D(T *A, int N)
{
	int ii;
	for (ii = 0; ii < N; ii++)
	{
		A[ii] /= N;
	}
}

class Limits_2D
{
	public:
		int XL, XU, YL, YU; //x and y, upper and lower

		//Constructor
		Limits_2D ()
		{
			XL = 0;
			YL = 0;
			XU = 1;
			YU = 1;
		}
};

//	Computational Parameters
class Comp_Parameters
{
	public:
		int procs_single,procs,Nx,Ny,nd,ne,eigen,normalize_wvfxn;
		double xmin,xmax,ymin,ymax,tolerance,dt,runtime;

		//Constructor
		Comp_Parameters ()
		{
			//Default values
			procs_single 		= 1;
			procs 				= 4;
			Nx 					= 8192;
			Ny 					= 1024;
			xmin 				= -400.0;
			xmax 				= 200.0;
			ymin 				= 9.0;
			ymax 				= 16.0;
			nd 					= 5;
			ne 					= 2;
			eigen 				= 1;
			tolerance 			= 0.005; //this might not be needed anymore
			normalize_wvfxn 	= 0;
			dt 					= 0.5;
			runtime			 	= 3000;
		}
};

//	Physical System Parameters
class Phys_Parameters
{
	public:
		double m,M,hbar,e_affinity,xi,kvect,Lbarrier,Rbarrier,work_fxn;
		double V_b,z0,epsilon,r_equil,rad,sigma,del,cconst,z_R,z_L;
		double omega,Lbarrier_itn,Rbarrier_itn,alpha,e_charge;
		int coupling;

		Phys_Parameters ()
		{
			//Default values
			m 			= 1.0;
			M 			= 2.59e6;
			hbar 		= 1.0;
			z_R 		= 24.57;
			z_L 		= 0.0;
			e_affinity 	= -2.56;
			xi 			= 40.0;
			kvect 		= 0.606;
			Lbarrier 	= 100.0;
			Rbarrier 	= 100.0;
			work_fxn 	= -5.4;
			V_b 		= 0.0;
			z0 			= -100.0;
			epsilon 	= 1.09;
			r_equil 	= 12.00;
			rad 		= 6.7085;
			sigma 		= 2.0;
			del 		= 6.2;
			coupling 	= 1;
			cconst 		= 1.0;
			alpha 		= 4.0;
			e_charge 	= -1.0;
		}
};

//	Output Data Parameters
class Out_Parameters
{
	public:
		int Potential,WvfxnInit,ElecStates,PhonStates,Animation,Slice,Correlation,FullWvfxn,Cnm;
		int SliceIndex,ElecCorr,MolecCorr,DCoeffR,DCoeffL,Phases,Energy,Lifetime,ZExp,FluxR,FluxL,FluxInt;

		Out_Parameters ()
		{
			//Default values
			Potential	=0;
			WvfxnInit	=0;
			ElecStates	=0;
			PhonStates	=0;
			Animation	=0;
			Slice		=0;
			SliceIndex	=0;
			ElecCorr	=0;
			MolecCorr	=0;
			DCoeffR		=0;
			DCoeffL		=0;
			Phases		=0;
			Energy		=0;
			Lifetime	=0;
			ZExp		=0;
			FluxR		=0;
			FluxL		=0;
			FluxInt		=0;
			Correlation	=0;
			FullWvfxn	=0;
			Cnm			=0;  
		}

};

//	Contains a series of 1D arrays to be used in ITP calculations
//	This class depends on functions contained in numerics.h
template <typename T> class EigenArray_1D
{
	public:
		int Nx; //number of grid points
		int n; //number of eigenstates
		T xstep;
		T *arrays;

		//Constructor
		EigenArray_1D (int xx, int yy)
		{
			Nx 		= xx;
			n 		= yy;
			arrays 	= new T [xx*yy];
		}

		//Destructor
		~EigenArray_1D ()
		{
			delete [] arrays;
		}

		//returns the memory address of the iith array in the set
		T *get_array_addr(int ii)
		{
			return arrays + ii*Nx;
		}

		//normalizes all arrays as wavefunctions, assuming phi*phi integrates to 1.00
		void normalize()
		{
			int ii,jj;
			T *num;
			Array_1D <T> A(Nx);
			A.xstep = real((cplx)xstep);
			
			for (ii = 0; ii < n; ii++)
			{
				num = get_array_addr(ii);
				for (jj = 0; jj < Nx; jj++)
				{
					A.grid[jj] = num[jj];
				}
				normalize_wxfxn_1D(A);
				for (jj = 0; jj < Nx; jj++)
				{
					num[jj] = A.grid[jj];
				}
			}
		}

		double normalize_n(int nn)
		{
			int ii;
			T *num;
			Array_1D <T> A(Nx);
			A.xstep = real(cplx(xstep));
			num 	= get_array_addr(nn);
			
			for (ii = 0; ii < Nx; ii++)
			{
				A.grid[ii] = num[ii];
			}
			double conv = real(integrate_1D(A));
			normalize_wxfxn_1D(A);
			for (ii = 0; ii < Nx; ii++)
			{
				num[ii] = A.grid[ii];
			}
			return conv;
		}
};
