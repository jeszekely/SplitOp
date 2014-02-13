#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <new>
#include <fftw3.h>
#include <string>
#include <vector>
#include <memory>
#include <boost/algorithm/string.hpp>
#include "iniparser.h"
#include "array_structs.h"
#include "numerics.h"
#include "junction_hamiltonian.h"
#include "imag_time_prop.h"
#include "quantum_dynamics.h"
#include <mkl.h>

using namespace std;

typedef std::complex<double> cplx;

int main () {

	int ii, jj, kk; //Variables useful for loops
	cplx comp;

/*************************************************************
	Import parameters from the input file (inputs.ini)
*************************************************************/

	Comp_Parameters Comp_Param;
	Phys_Parameters Phys_Param;
	Out_Parameters Out_Param;

 	dictionary * ini = iniparser_load("inputs.ini");
    if (ini==NULL) {
        fprintf(stderr, "Cannot parse file: inputs.ini\n");
        return -1 ;
    }
    iniparser_dump(ini, stderr);

	Comp_Param.procs 			= iniparser_getint(ini,"CompCell:Procs",4);
    Comp_Param.xmin 			= iniparser_getdouble(ini, "CompCell:Xmin", -400.0);
    Comp_Param.xmax 			= iniparser_getdouble(ini, "CompCell:Xmax", 200.0);
	Comp_Param.ymin 			= iniparser_getdouble(ini, "CompCell:Ymin", 9.0);
    Comp_Param.ymax 			= iniparser_getdouble(ini, "CompCell:Ymax", 16.0);
    Comp_Param.Nx 				= iniparser_getint(ini, "CompCell:Nx", 8192);
    Comp_Param.Ny 				= iniparser_getint(ini, "CompCell:Ny", 1064);
    Comp_Param.nd 				= iniparser_getint(ini,"CompCell:Nd", 5);
    Comp_Param.ne 				= iniparser_getint(ini,"CompCell:Ne", 3);
    Comp_Param.eigen 			= iniparser_getint(ini,"CompCell:Eigen",1);
    Phys_Param.z_R 				= iniparser_getdouble(ini,"CompCell:zR",24.57);
	Phys_Param.z_L 				= iniparser_getdouble(ini,"CompCell:zL",0.0);
    Comp_Param.tolerance 		= iniparser_getdouble(ini, "CompCell:Tolerance", 0.01);
    Comp_Param.normalize_wvfxn 	= iniparser_getint(ini,"CompCell:Normalize",0);

    Comp_Param.dt 				= iniparser_getdouble(ini, "TimeParameters:Dt", 0.5);
    Comp_Param.runtime 			= iniparser_getdouble(ini, "TimeParameters:Runtime",3000.0);

	Phys_Param.hbar 		= iniparser_getdouble(ini,"PhysicalParameters:hbar",1.0);
    Phys_Param.m 			= iniparser_getdouble(ini, "PhysicalParameters:M_electron", 1.0);
	Phys_Param.M 			= iniparser_getdouble(ini, "PhysicalParameters:M_C60", 2.59e6);
	Phys_Param.e_affinity 	= iniparser_getdouble(ini,"PhysicalParameters:EAffin",-2.65);
	Phys_Param.xi 			= iniparser_getdouble(ini, "PhysicalParameters:Xi", 40.0);
	Phys_Param.kvect 		= iniparser_getdouble(ini, "PhysicalParameters:ElecVel", 0.606);
	Phys_Param.Lbarrier 	= iniparser_getdouble(ini, "PhysicalParameters:LBarrier", 100.0);
	Phys_Param.Rbarrier 	= iniparser_getdouble(ini, "PhysicalParameters:RBarrier", 100.0);
	Phys_Param.work_fxn 	= iniparser_getdouble(ini, "PhysicalParameters:WorkFxn", -5.4);
	Phys_Param.V_b 			= iniparser_getdouble(ini, "PhysicalParameters:Bias", 0.0);
	Phys_Param.z0 			= iniparser_getdouble(ini, "PhysicalParameters:ElecPos", -150.0);
	Phys_Param.epsilon 		= iniparser_getdouble(ini, "PhysicalParameters:Epsilon", 1.09);
	Phys_Param.r_equil 		= iniparser_getdouble(ini, "PhysicalParameters:Requil", 11.57);
	Phys_Param.rad 			= iniparser_getdouble(ini, "PhysicalParameters:Rad", 6.7085);
	Phys_Param.sigma 		= iniparser_getdouble(ini, "PhysicalParameters:Sigma",  2.0);
	Phys_Param.del 			= iniparser_getdouble(ini, "PhysicalParameters:Del", 6.2);
	Phys_Param.coupling 	= iniparser_getint(ini, "PhysicalParameters:CoupFxn", 1);
	Phys_Param.cconst 		= iniparser_getdouble(ini,"PhysicalParameters:CConst",1.0);


	Out_Param.Potential 	= iniparser_getint(ini, "Outputs:Potential", 0);
	Out_Param.WvfxnInit 	= iniparser_getint(ini, "Outputs:WvfxnInit", 0);
	Out_Param.ElecStates 	= iniparser_getint(ini, "Outputs:ElecStates", 0);
	Out_Param.PhonStates 	= iniparser_getint(ini, "Outputs:PhonStates", 0);
	Out_Param.Animation 	= iniparser_getint(ini, "Outputs:Animation", 0);
	Out_Param.Slice 		= iniparser_getint(ini, "Outputs:Slice", 0);
	Out_Param.SliceIndex 	= iniparser_getint(ini, "Outputs:SliceIndex", 0);
	Out_Param.ElecCorr 		= iniparser_getint(ini, "Outputs:ElecCorr", 0);
	Out_Param.MolecCorr 	= iniparser_getint(ini, "Outputs:MolecCorr", 0);
	Out_Param.DCoeffR 		= iniparser_getint(ini, "Outputs:DCoeffR", 0);
	Out_Param.DCoeffL 		= iniparser_getint(ini, "Outputs:DCoeffL", 0);
	Out_Param.Phases 		= iniparser_getint(ini, "Outputs:Phases", 0);
	Out_Param.Energy 		= iniparser_getint(ini, "Outputs:Energy", 0);
	Out_Param.Lifetime 		= iniparser_getint(ini, "Outputs:Lifetime", 0);
	Out_Param.ZExp 			= iniparser_getint(ini, "Outputs:ZExp", 0);
	Out_Param.FluxR			= iniparser_getint(ini, "Outputs:FluxR", 0);
	Out_Param.FluxL 		= iniparser_getint(ini, "Outputs:FluxL", 0);
	Out_Param.FluxInt 		= iniparser_getint(ini, "Outputs:FluxInt", 0);
	Out_Param.Correlation 	= iniparser_getint(ini, "Outputs:Correlation", 0);
	Out_Param.FullWvfxn 	= iniparser_getint(ini, "Outputs:FullWvfxn", 0);
	Out_Param.Cnm 			= iniparser_getint(ini, "Outputs:Cnm",0); 

    iniparser_freedict(ini);

/*************************************************************
	Calculate derived parameters
*************************************************************/

    Phys_Param.epsilon 		/= 27.21; //Converts epsilon into atomic units
    Phys_Param.e_affinity 	/= 27.21; //Converts electron affinity into atomic units
    Phys_Param.work_fxn 	/= 27.21; //Converts work function into atomic units

    Phys_Param.omega 		= sqrt(72 * Phys_Param.epsilon/(Phys_Param.r_equil*Phys_Param.r_equil) * (1/Phys_Param.M)); //in au
    Phys_Param.Lbarrier_itn = 2*M_PI*M_PI/(Phys_Param.Lbarrier*Phys_Param.Lbarrier);
	Phys_Param.Rbarrier_itn = 2*M_PI*M_PI/(Phys_Param.Rbarrier*Phys_Param.Rbarrier);

/*************************************************************
	Initialize the fftw threading
*************************************************************/

	fftw_init_threads(); //initialize SMP
	fftw_plan_with_nthreads(Comp_Param.procs); //use #(procs) processors

/*************************************************************
	Define Potential Arrays for the main 2D propagation
*************************************************************/

//	Spatial grid step sizes
	double xstep = (Comp_Param.xmax - Comp_Param.xmin)/Comp_Param.Nx;
	double ystep = (Comp_Param.ymax - Comp_Param.ymin)/Comp_Param.Ny;
	double pstep = M_PI * 2.0 /(Comp_Param.Nx * xstep);
	double qstep = M_PI * 2.0 /(Comp_Param.Ny * ystep);

	cout << "The grid step sizes are:" << endl << "xstep = " << xstep << endl << "ystep = " << ystep << endl << "pstep = " << pstep << endl <<"qstep = " << qstep << endl;

//	Define 1D position and momentum space arrays
	Array_1D <double> Xgrid(Comp_Param.Nx);
	Xgrid.xinit = Comp_Param.xmin;
	Xgrid.xstep = xstep;
	Xgrid.fill_array();

	Array_1D <double> Ygrid(Comp_Param.Ny);
	Ygrid.xinit = Comp_Param.ymin;
	Ygrid.xstep = ystep;
	Ygrid.fill_array();

	Array_1D <double> Pgrid(Comp_Param.Nx);
	Pgrid.xinit = 0.0;
	Pgrid.xstep = pstep;
	Pgrid.fill_array();

	Array_1D <double> Qgrid(Comp_Param.Ny);
	Qgrid.xinit = 0.0;
	Qgrid.xstep = qstep;
	Qgrid.fill_array();

//	Correct the ordering of the momentum space arrays for the fftw algorithm
	for (ii = Comp_Param.Nx/2; ii < Comp_Param.Nx; ii++)
	{
		Pgrid.grid[ii] = -1.0*Pgrid.grid[Comp_Param.Nx - ii];
	}

	for (jj = Comp_Param.Ny/2; jj < Comp_Param.Ny; jj++)
	{
		Qgrid.grid[jj] = -1.0*Qgrid.grid[Comp_Param.Ny - jj];
	}

//	2D Potential and Kinetic energy arrays
	Array_2D <cplx> V(Comp_Param.Nx,Comp_Param.Ny);
	V.xstep = xstep;
	V.ystep = ystep;

	for (ii = 0; ii < V.Nx; ii++)
	{
		for (jj = 0; jj < V.Ny; jj++)
		{
			V.set_elem(ii,jj,H_Potential(Xgrid.grid[ii],Ygrid.grid[jj],Phys_Param,Comp_Param));
		}
	}

	Array_2D <cplx> T(Comp_Param.Nx,Comp_Param.Ny);
	T.xstep = pstep;
	T.ystep = qstep;

	for (ii = 0; ii < T.Nx; ii++)
	{
		for (jj = 0; jj < T.Ny; jj++)
		{
			T.set_elem(ii,jj,H_Kinetic(Pgrid.grid[ii],Qgrid.grid[jj],Phys_Param));
		}
	}

	Array_2D <cplx> KinetOp(T.Nx,T.Ny);
	for (ii=0; ii < T.Nx*T.Ny; ii++)
	{
		comp				= -0.5j * Comp_Param.dt * T.grid[ii];
		KinetOp.grid[ii] 	= exp(comp);
	}

	Array_2D <cplx> PotenOp(V.Nx,V.Ny);
	for (ii=0; ii < V.Nx*V.Ny; ii++)
	{
		comp 				= -1.0j * Comp_Param.dt * V.grid[ii];
		PotenOp.grid[ii] 	= exp(comp);
	}

/*************************************************************
	Define Potential Arrays for the 1D Calculations
*************************************************************/

//	1D Molecular Arrays for ITP
	Array_1D <double> V_nuc(Comp_Param.Ny);
	V_nuc.xstep = ystep;

	Array_1D <double> T_nuc(Comp_Param.Ny);
	T_nuc.xstep = qstep;

	for (ii = 0; ii < Comp_Param.Ny; ii++)
	{
		V_nuc.grid[ii] = V_M(Ygrid.grid[ii],Phys_Param);
		T_nuc.grid[ii] = H_Kinetic(0.0,Qgrid.grid[ii],Phys_Param);
	}

//	1D Electronic Arrays for ITP
	int JL = -1; //JL an JR are the array indexes of Xgrid the correspond to the junction region for ITP
	int JR = -1;
	for (ii = 0; ii < Comp_Param.Nx; ii++)
	{
		if ((JL < 0) && (Xgrid.grid[ii] >= (Phys_Param.r_equil - 45.0) ))
		{
			JL = ii;
		}
		if ((JR < 0) && (Xgrid.grid[ii] >= (Phys_Param.r_equil + 45.0) ))
		{
			JR = ii;
		}
		if (JL > 0 && JR > 0)
		{
			break;
		}
	}

	int Nx_junc = JR - JL;
	if ((Nx_junc % 2) == 1) //Adjusts the region to contain an even number of grid points if necessary
	{
		Nx_junc++;
		JL++;
	}

	cout << "The junction spans from index " << JL << " to index " << JR << " in the electronic dimension." << endl;
	cout << flush;

//	Position grid in the junction region
	Array_1D <double> Xgrid_junc(Nx_junc);
	Xgrid_junc.xstep = xstep;
	Xgrid_junc.xinit = Xgrid.grid[JL];
	Xgrid_junc.fill_array();

//	Momentum grid in the junction region
	Array_1D <double> Pgrid_junc(Nx_junc);
	double pstep_junc = M_PI * 2.0 / (xstep*Nx_junc);
	for (ii = 0; ii < Nx_junc; ii++)
	{
		Pgrid_junc.grid[ii] = ii*pstep_junc;
	}
	for (ii = Nx_junc/2; ii < Nx_junc; ii++)
	{
		Pgrid_junc.grid[ii] = -1.0*Pgrid_junc.grid[Nx_junc - ii];
	}

//	Energy Arrays
	Array_1D <double> PJE (Nx_junc);
	PJE.xstep = xstep;
	Array_1D <double> KJE (Nx_junc);
	KJE.xstep = pstep_junc;

	for (ii = 0; ii < Nx_junc; ii++)
	{
		PJE.grid[ii] = Elec_parab(Xgrid_junc.grid[ii],Phys_Param.r_equil,Phys_Param);
		KJE.grid[ii] = Pgrid_junc.grid[ii] * Pgrid_junc.grid[ii] / (2.0 * Phys_Param.m);
	}

/*************************************************************
	Imaginary Time Propagation in the Electronic Subspace
*************************************************************/

	cout << "Getting Electron Eigenstates" << endl;

	EigenArray_1D <cplx> ElecStates(Nx_junc,Comp_Param.ne);
	ElecStates.xstep = xstep;
	get_wvfxns(ElecStates,Xgrid_junc,PJE,KJE,Phys_Param,Phys_Param.m);

//	Print results to file
	if (Out_Param.ElecStates == 1)
	{
		ofstream ElecITP;
		ElecITP.open("output_data/ElecITP_Results.txt");
		for (kk = 0; kk < Nx_junc; kk++)
		{
			ElecITP << Xgrid_junc.grid[kk] << "\t" << PJE.grid[kk];
			for (jj = 0; jj < ElecStates.n; jj++)
			{
				ElecITP << "\t" << real((ElecStates.get_array_addr(jj))[kk]);
			}
			ElecITP << endl;
		}
		ElecITP.close();
	}


/*************************************************************
	Imaginary Time Propagation in the Molecular Subspace
*************************************************************/

	cout << "Getting Phonon Eigenstates" << endl;
	EigenArray_1D <cplx> MolecStates(Comp_Param.Ny,Comp_Param.nd);
	MolecStates.xstep = ystep;
	if (Comp_Param.eigen == 1)
	{
		get_wvfxns(MolecStates,Ygrid,V_nuc,T_nuc,Phys_Param,Phys_Param.M);
	}
	else if (Comp_Param.eigen == 2)
	{
		load_eigen_array(MolecStates);
	}
	else
	{
		cplx *wvfxn_ptr;
		for (ii = 0; ii < Comp_Param.nd; ii++)
		{
			wvfxn_ptr = MolecStates.get_array_addr(ii);
			for (jj = 0; jj < MolecStates.Nx; jj++)
			{
				wvfxn_ptr[jj] = (cplx)ho_wvfxn(ii,Ygrid.grid[jj],Phys_Param.r_equil, Phys_Param.omega, Phys_Param.M, Phys_Param.hbar);
			}
		}
	}

	if (Out_Param.PhonStates == 1)
	{
		ofstream MolITP;
		MolITP.open("output_data/MolITP_Results.txt");
		for (kk = 0; kk < MolecStates.Nx; kk++)
		{
			MolITP << Ygrid.grid[kk]<< "\t" << V_nuc.grid[kk];
			for (jj = 0; jj < MolecStates.n; jj++)
			{
				MolITP << "\t" << real((MolecStates.get_array_addr(jj))[kk]);
			}
			MolITP << endl;
		}
		MolITP.close();
	}

/***********************************************************************
	Determine where the imaginary potential begins and ends (electronic)
************************************************************************/

	int Poten_L = 0;
	int Poten_R = 0;
	ii = 0;
	int found = 0;
	while (found == 0)
	{
		if (imag(V.get_elem(Poten_L,Comp_Param.Ny/2)) == 0)
		{
			found = 1;
		}
		else
		{
			Poten_L++;
		}
	}
	cout << "Left electronic absorbing potential ends at array index " << Poten_L << "." << endl;

	if (Poten_L > Xgrid.Nx)
	{
		cout << "I've made a huge mistake." << endl;
		exit(-1);
	}
	Poten_R = Poten_L;
	while (found == 1)
	{
		if (imag(V.get_elem(Poten_R,Comp_Param.Ny/2)) != 0)
		{
			found = 0;
		}
		else
		{
			Poten_R++;
		}
	}
	cout << "Right electronic absorbing potential begins at array index " << Poten_R << "." << endl;
	if (Poten_R > Xgrid.Nx || Poten_R < Poten_L)
	{
		cout << "I've made a huge mistake." << endl;
		exit(-1);
	}

//	Locations to calculate flux, moved away from imaginary surface for calculation of numerical derivatives
	int zR = Poten_R - 10;
	int zL = Poten_L + 10;

/***********************************************************************
	Determine where the imaginary potential begins and ends (electronic)
************************************************************************/

	int MolPoten_L = 0;
	int MolPoten_R = 0;
	ii = 0;
	found = 0;
	while (found == 0)
	{
		if (imag(V.get_elem(zR,MolPoten_L)) == 0)
		{
			found = 1;
		}
		else
		{
			MolPoten_L++;
		}
	}
	cout << "Left molecular absorbing potential ends at array index " << MolPoten_L << "." << endl;

	if (MolPoten_L > Xgrid.Nx)
	{
		cout << "I've made a huge mistake." << endl;
		exit(-1);
	}
	MolPoten_R = MolPoten_L;
	while (found == 1)
	{
		if (imag(V.get_elem(zR,MolPoten_R)) != 0)
		{
			found = 0;
		}
		else
		{
			MolPoten_R++;
		}
	}
	cout << "Right molecular absorbing potential begins at array index " << MolPoten_R << "." << endl;
	if (MolPoten_R > Ygrid.Nx || MolPoten_R < MolPoten_L)
	{
		cout << "I've made a huge mistake." << endl;
		exit(-1);
	}

//	Locations to calculate flux, moved away from imaginary surface for calculation of numerical derivatives
	int MolzR = MolPoten_R - 5;
	int MolzL = MolPoten_L + 5;


/*************************************************************
	Determine where the junction begins and ends
*************************************************************/

	int zR_junc, zL_junc;
	for (ii = 0; ii < Xgrid.Nx; ii++)
	{
		if (Xgrid.grid[ii] >= Phys_Param.z_L)
		{
			zL_junc = ii;
			break;
		}
	}
	for (ii = 0; ii<Xgrid.Nx; ii++)
	{
		if (Xgrid.grid[ii] >= Phys_Param.z_R)
		{
			zR_junc = ii;
			break;
		}
	}
	cout << "The junction spans from index " << zL_junc << " to index " << zR_junc << "." << endl;

/*************************************************************
	Define Wavefunction Array, Fourier Transform setup
*************************************************************/

	Array_2D <cplx> Wvfxn(Comp_Param.Nx,Comp_Param.Ny);

//	Forward and backward Fourier transforms
	fftw_plan_with_nthreads(Comp_Param.procs); //reset since ITP may have altered this value
    fftw_plan forplan, backplan;

	forplan = fftw_plan_dft_2d(Wvfxn.Nx,Wvfxn.Ny,(fftw_complex *)Wvfxn.grid,\
	(fftw_complex*)Wvfxn.grid, FFTW_FORWARD, FFTW_MEASURE);

	backplan = fftw_plan_dft_2d(Wvfxn.Nx,Wvfxn.Ny,(fftw_complex *)Wvfxn.grid,\
	(fftw_complex*)Wvfxn.grid, FFTW_BACKWARD, FFTW_MEASURE);

	Wvfxn.xstep = xstep;
	Wvfxn.ystep = ystep;

	cplx *wvfxn_nuc;
	wvfxn_nuc = MolecStates.get_array_addr(0); //Uses the ground state phonon mode
	for (ii = 0; ii < Wvfxn.Nx; ii++)
	{
		for (jj = 0; jj < Wvfxn.Ny; jj++)
		{
			Wvfxn.set_elem(ii,jj,wvfxn_electron(Xgrid.grid[ii],Phys_Param)*wvfxn_nuc[jj]);
		}
	}
	normalize_wxfxn_2D(Wvfxn);

//	Make a copy array for the correlation funciton calculation
	Array_2D <cplx> Wvfxn_Init(Comp_Param.Nx,Comp_Param.Ny);
	Wvfxn_Init.xstep = xstep;
	Wvfxn_Init.ystep = ystep;
	for (ii = 0; ii < Comp_Param.Nx*Comp_Param.Ny; ii++)
	{
		Wvfxn_Init.grid[ii] = Wvfxn.grid[ii];
	}

//	Allocates array space for use in several subroutines that would
//	otherwise have to allocate memory each time they are called

	Array_2D <cplx> Scratch(Comp_Param.Nx,Comp_Param.Ny);
	Scratch.xstep = xstep;
	Scratch.ystep = ystep;

/*************************************************************
	Output Any Initial conditions described in inputs.ini
*************************************************************/

//	Output the full 2D potential energy surface
	if (Out_Param.Potential == 1)
	{
		cout << "Outputting the full 2D Potential energy surface..."<< endl;
		ofstream PotenFile;
		PotenFile.open("output_data/Potential2D.txt");

		for (ii = 0; ii < V.Nx; ii += 10) //step by 10 otherwise output file is enormous
		{
			for (jj = 0; jj < V.Ny; jj += 4)
			{
				PotenFile << Xgrid.grid[ii] << "\t";
				PotenFile << Ygrid.grid[jj] << "\t";
				PotenFile << real(V.get_elem(ii,jj)) << "\t" << imag(V.get_elem(ii,jj)) << endl;
			}
			PotenFile << endl;
		}
		PotenFile.close();
	}

//	Output the full 2D initial wavefunction
	if (Out_Param.WvfxnInit == 1)
	{
		cout << "Outputting the full 2D Wavefunction..." << endl;
		ofstream WvfxnInitFile;
		WvfxnInitFile.open("output_data/Wavefunction_Initial.txt");
		for (ii = 0; ii < Wvfxn.Nx; ii += 10) //step by 10 otherwise output file is enormous
		{
			for (jj = 0; jj < Wvfxn.Ny; jj += 4)
			{
				WvfxnInitFile << Xgrid.grid[ii] << "\t";
				WvfxnInitFile << Ygrid.grid[jj] << "\t";
				WvfxnInitFile << real(Wvfxn_Init.get_elem(ii,jj)) << "\t" << imag(Wvfxn_Init.get_elem(ii,jj));
				WvfxnInitFile << scientific << endl;
			}
			WvfxnInitFile << endl;
		}

		WvfxnInitFile.close();
	}

/*************************************************************
	Allocate small arrays for the various calculations
*************************************************************/

//	Arrays in which the dcoefficient values are returned
	Array_1D <cplx> DCoeff_L(Comp_Param.nd);
	Array_1D <cplx> DCoeff_R(Comp_Param.nd);

//	Array to return the values from the energy calculation
	Array_1D <double> Energy(5); //contains ElecKE, ElecPE, MolKE, MolPE and the interaction
	Array_1D <double> JuncEnergy(5); //Energy components just in the junction region

//	Limits for the Junction Wavefunction Calculation
	Limits_2D Junction_Limits;
	Junction_Limits.XU = zR_junc;
	Junction_Limits.XL = zL_junc;
	Junction_Limits.YU = MolPoten_R - 1;
	Junction_Limits.YL = MolPoten_L + 1;

//	Limits for normal integration, ensures imaginary potential is not included in the energy integration
	Limits_2D CompCellLimits;
	CompCellLimits.XU = Poten_R - 1;
	CompCellLimits.XL = Poten_L + 1;
	CompCellLimits.YU = MolPoten_R - 1;
	CompCellLimits.YL = MolPoten_L + 1;

/****************************************************************
	Open the necessary output files
****************************************************************/

	ofstream EnergyFile;
	ofstream JuncEnergyFile;
	ofstream LifetimeFile;
	ofstream ZExpFile;
	ofstream CorrelationFile;
	ofstream DCoeffRFile;
	ofstream DCoeffLFile;
	ofstream FluxRFile;
	ofstream FluxLFile;
	ofstream WvfxnOut;
	ofstream CnmOut; 

	if (Out_Param.Energy != 0)
	{
		EnergyFile.open("output_data/Energy.txt");
		EnergyFile << "#Time ElectronKE MoleculeKE ElectronPE MoleculePE InteractionEnergy" << endl;

		JuncEnergyFile.open("output_data/JunctionEnergy.txt");
		JuncEnergyFile << "#Time ElectronKE MoleculeKE ElectronPE MoleculePE InteractionEnergy" << endl;
	}
	if (Out_Param.Lifetime != 0)
	{
		LifetimeFile.open("output_data/Lifetime.txt");
		LifetimeFile << "#Time\t Probability Density in Junction" << endl;
	}
	if (Out_Param.ZExp != 0)
	{
		ZExpFile.open("output_data/ZExp.txt");
		ZExpFile << "#Time\t Wvfxn Norm\t Molecule(center of mass) expectation value" << endl;
	}
	if (Out_Param.Correlation != 0)
	{
		CorrelationFile.open("output_data/Correlation.txt");
		CorrelationFile << "#Time\t Correlation" << endl;
	}
	if (Out_Param.DCoeffR != 0)
	{
		DCoeffRFile.open("output_data/DCoeff_Transmitted.txt");
		DCoeffRFile << "#Time";
		for (ii = 0; ii < DCoeff_R.Nx; ii++)
		{
			DCoeffRFile << "\t" << ii;
		}
		DCoeffRFile << endl;
	}
	if (Out_Param.DCoeffL != 0)
	{
		DCoeffLFile.open("output_data/DCoeff_Reflected.txt");
		DCoeffLFile << "#Time";
		for (ii = 0; ii < DCoeff_L.Nx; ii++)
		{
			DCoeffLFile << "\t" << ii;
		}
		DCoeffLFile << endl;
	}
	if (Out_Param.FluxR != 0)
	{
		FluxRFile.open("output_data/Flux_Transmitted.txt");
		FluxRFile << "#Time\t Transmitted Flux" << endl;
	}
	if (Out_Param.FluxL != 0)
	{
		FluxLFile.open("output_data/Flux_Reflected.txt");
		FluxLFile << "#Time\t Reflected Flux" << endl;
	}

	if (Out_Param.Cnm != 0)
	{
		CnmOut.open("output_data/Cnm.txt");
		CnmOut << "#"; 
		for (ii = 0; ii < ElecStates.n; ii++)
		{
			for (jj = 0; jj < MolecStates.n; jj++)
			{
				CnmOut << "C" << ii << jj << " ";  
			}
		}
		CnmOut << endl; 
	}
/****************************************************************
	Primary calculation loop, steps a wavefunction forward in time
****************************************************************/

	char filename[sizeof("output_data/2Dwvfxn_t000000.txt")];
	int tindex = 0;
	cplx corr;
	do
	{
//	Output full 2D wavefunction
		cout << "Step #" << tindex << " time = " << scientific << Wvfxn.time << endl;
		if (Out_Param.Animation != 0)
		{
			if (tindex % Out_Param.Animation == 0)
			{
				print_wvfxn(Wvfxn,Xgrid,Ygrid,tindex);
			}
		}

//	Output full wavefunction during propagation
	if (Out_Param.FullWvfxn != 0)
	{
		if (tindex % Out_Param.FullWvfxn == 0)
		{
			sprintf(filename,"output_data/2Dwvfxn_t%.6d.txt",tindex); //Appended filename
		    WvfxnOut.open(filename);
	   		for (ii = 0; ii < Wvfxn.Nx; ii += 10)
   			{
				for (jj = 0; jj < Wvfxn.Ny; jj += 4)
				{
					WvfxnOut << Xgrid.grid[ii] << "\t" << Ygrid.grid[jj] << "\t" << abs(Wvfxn.get_elem(ii,jj)) << endl;
				}
				WvfxnOut << endl;
			}
			WvfxnOut.close();

		}
	}

//	Output wavefunction slices
		if (Out_Param.Slice != 0)
		{
			if (tindex % Out_Param.Slice == 0)
			{
				print_wvfxn_slice(Wvfxn,Xgrid,Out_Param.SliceIndex,tindex);
			}
		}

//	Calculate the energy parameters of the system, print them to file
		if (Out_Param.Energy != 0)
		{
			if (tindex % Out_Param.Energy == 0)
			{
				get_energy_components(Wvfxn, Xgrid, Ygrid, Energy, JuncEnergy, CompCellLimits, Junction_Limits, Phys_Param, Scratch);

				EnergyFile 	<< Wvfxn.time << "\t" << Energy.grid[0] << "\t" << Energy.grid[1]\
							<< "\t" << Energy.grid[2] << "\t" << Energy.grid[3] << "\t" << Energy.grid[4] << endl;

				JuncEnergyFile 	<< Wvfxn.time << "\t" << JuncEnergy.grid[0] << "\t" << JuncEnergy.grid[1]\
								<< "\t" << JuncEnergy.grid[2] << "\t" << JuncEnergy.grid[3] << "\t" << JuncEnergy.grid[4] << endl;
			}
		}

//	Calculate the wavefunction density in the junction
		if (Out_Param.Lifetime != 0)
		{
			if (tindex % Out_Param.Lifetime == 0)
			{
				LifetimeFile 	<< Wvfxn.time << "\t"\
								<< get_junction_wavefunction(Wvfxn,Junction_Limits,Scratch) << endl;
			}
		}

//	Calculate the the molecular expectation value
		if (Out_Param.ZExp != 0)
		{
			if (tindex % Out_Param.ZExp == 0)
			{
				ZExpFile << Wvfxn.time << "\t" << get_wvfxn_norm_2D(Wvfxn,CompCellLimits,Scratch) \
				<< "\t" << Z_ExpVal_2D(Wvfxn,Ygrid,CompCellLimits,Scratch) << endl;
			}
		}

//	Calculate the Correlation function
		if (Out_Param.Correlation != 0)
		{
			if (tindex % Out_Param.Correlation == 0)
			{
				corr = correlation(Wvfxn,Wvfxn_Init,Scratch);
				CorrelationFile << Wvfxn.time << "\t" << corr.real() << "\t" << corr.imag() << endl;
			}
		}

//	Calculate the transmitted DCoefficient
		if (Out_Param.DCoeffR != 0)
		{
			get_Dcoefficients(Wvfxn,MolecStates,DCoeff_R,zR);
			DCoeffRFile << Wvfxn.time;
			for (ii = 0; ii < DCoeff_R.Nx; ii++)
			{
				DCoeffRFile << "\t" << DCoeff_R.grid[ii].real() << "\t" << DCoeff_R.grid[ii].imag();
			}
			DCoeffRFile << endl;
		}

//	Calculate the reflected DCoefficient
		if (Out_Param.DCoeffL != 0)
		{
			get_Dcoefficients(Wvfxn,MolecStates,DCoeff_L,zL);
			DCoeffLFile << Wvfxn.time;
			for (ii = 0; ii < DCoeff_L.Nx; ii++)
			{
				DCoeffLFile << "\t" << DCoeff_L.grid[ii].real() << "\t" << DCoeff_L.grid[ii].imag();
			}
			DCoeffLFile << endl;
		}

//	Calculate the transmitted flux
		if (Out_Param.FluxR != 0)
		{
			if (tindex % Out_Param.FluxR == 0)
			{
				FluxRFile << Wvfxn.time << "\t" << probcurrent_2D(Wvfxn,zR,Phys_Param.m, Phys_Param.hbar) << endl;
			}
		}

//	Calculate the reflected flux
		if (Out_Param.FluxL != 0)
		{
			if (tindex % Out_Param.FluxL == 0)
			{
				FluxLFile << Wvfxn.time << "\t" << -1.0*probcurrent_2D(Wvfxn,zL,Phys_Param.m, Phys_Param.hbar) << endl;
			}
		}

		if (Out_Param.Cnm != 0)
		{
			if (tindex % Out_Param.Cnm == 0)
			{
				CnmOut << Wvfxn.time; 
				for (ii = 0; ii < ElecStates.n; ii++)
				{
					for (jj = 0; jj < MolecStates.n; jj++)
					{
						CnmOut << " " << abs(get_junction_Cnm(Wvfxn,ElecStates,MolecStates,Junction_Limits,jj,ii,JL,Scratch));
					}
				}
				CnmOut << endl; 
			}
		}
//	Step the wavefunction and time parameters
		SplitOp2D_Step(Wvfxn,KinetOp,PotenOp,forplan,backplan);
		Wvfxn.time += Comp_Param.dt;
		tindex++;
	} while (Wvfxn.time < Comp_Param.runtime);

/*************************************************************
	Close output files
*************************************************************/

	if (Out_Param.Energy != 0) {EnergyFile.close(); JuncEnergyFile.close();}
	if (Out_Param.Lifetime != 0) {LifetimeFile.close();}
	if (Out_Param.ZExp != 0) {ZExpFile.close();}
	if (Out_Param.Correlation != 0) {CorrelationFile.close();}
	if (Out_Param.DCoeffR != 0) {DCoeffRFile.close();}
	if (Out_Param.DCoeffL != 0) {DCoeffLFile.close();}
	if (Out_Param.FluxR != 0) {FluxRFile.close();}
	if (Out_Param.FluxL != 0) {FluxLFile.close();}
	if (Out_Param.Cnm != 0)	{CnmOut.close();}


/*******************************************************************
	Below are functions to test the code while being written
	These are temporary and will be deleted when the code is complete
*********************************************************************/
/*
	cout << "Test #1\n";

	double *gridparams;

	gridparams = new double [4];
	gridparams[0] = 1.0;
	gridparams[1] = 1.0;
	gridparams[2] = 1.0;
	gridparams[3] = 1.0;

	Array_2D <cplx> array (10,10,gridparams); //Declare test array

	//Check get and set elements
	array.set_elem(5,5,3+2.0j);
	cout << "Test #2... the element at (5,5) is " << array.get_elem(5,5) << "\n";

	//Check addresses in memory
	cout << "Test #3... \n" <<  &array << "\n" << array.grid << "\n" << &array.grid << "\n" << &array.grid[1] << "\n";

	cout << sizeof(cplx) << "\n";

	delete [] gridparams;

	fftw_normalize_2D(array);
	cout << "Test #4... the element at (5,5) is now " << array.get_elem(5,5) << "\n";


	Array_1D<double> sine (200); //Declare test array
	for (ii = 0; ii < sine.Nx; ii++){
		sine.grid[ii] = sin(2*M_PI * ii / sine.Nx);
	}
	sine.xstep = (2*M_PI/sine.Nx);
	double x = deriv_1D_1st(sine,sine.Nx/2);
	double y = deriv_1D_2nd(sine,sine.Nx/2);
	cout.precision(10);
	cout << "Test #5... numeric derivative " << x << " " << y << endl;

	double integral = integrate_1D(sine);
	cout << "Test #6... numeric integration " << integral << endl;

	Array_1D <double> xgrid (200);
	Array_1D <double> ygrid (300);
	Array_2D <double> fxn (xgrid.Nx,ygrid.Nx);

	for (ii = 0; ii < xgrid.Nx; ii++){
		xgrid.grid[ii] = ii*2*M_PI/xgrid.Nx;
	}
	for (ii = 0; ii < ygrid.Nx; ii++){
		ygrid.grid[ii] = ii*2*M_PI/ygrid.Nx;
	}
	for (ii = 0; ii < xgrid.Nx; ii++) {
		for (jj = 0; jj < ygrid.Nx; jj++) {
			fxn.set_elem(ii,jj,sin(xgrid.grid[ii] + 2*ygrid.grid[jj]));
		}
	}
	fxn.xstep = xgrid.grid[1] - xgrid.grid[0];
	fxn.ystep = ygrid.grid[1] - ygrid.grid[0];

	double dx = partial_x(fxn,30,120);
	double dy = partial_y(fxn,30,120);
	double dx2 = cos(xgrid.grid[30] + 2*ygrid.grid[120]);
	double dy2 = 2*cos(xgrid.grid[30] + 2*ygrid.grid[120]);

	cout << "Test #7... 2D numeric derivatives " << endl << dx << endl << dx2 << endl << dy << endl << dy2 << endl;

	Limits_2D intlims;
	intlims.XL = 0;
	intlims.XU = xgrid.Nx/3;
	intlims.YL = 0;
	intlims.YU = ygrid.Nx/6;
	double integral2 = integrate_2D(fxn,intlims);

	cout << "Test #8... 2D numeric integration " << endl << integral2 << endl;

	cplx compnum = 1.0j* 0.25 * M_PI;
	compnum = exp (compnum);
	cout << "Test #9... complex exponential, e^(I*Pi/4) = " << compnum << endl;

	EigenArray_1D <double> eigen (Comp_Param.Nx,3);
	eigen.xstep = Xgrid.xstep;
	cout << "Test #10... testing EigenArray class" << endl;
	cout << "object address: " << &eigen << endl;
	cout << "object grid address: " << eigen.arrays << endl;
	cout << "for reference, a double has size " << sizeof(double) << endl;
	for (ii = 0; ii<3; ii++){
		cout << eigen.get_array_addr(ii) << endl;
	}

	cout << "Test #11... function initialization procedure" << endl;
	initialize_fxn_array(eigen,Xgrid);

	cplx proj = projection(eigen,1,2);
	cout << "The projection is: " << real(proj) << endl;
	eigen.normalize();
	proj = projection(eigen,1,2);
	cout << "The projection is now: " << real(proj) << " after normalization" << endl;

	load_eigen_array(eigen);

	cout << "Test #12... Imaginary time propagation check" << endl;
	Array_1D <double> Potential_1D_Test(2000);
	Array_1D <double> Kinetic_1D_Test(2000);

	double xstep_test = (10.0 - (-10.0))/2000;
	double pstep_test = M_PI * 2.0 /(2000 * xstep_test);

	Array_1D <double> Xgrid_test(2000);
	Xgrid_test.xinit = -10.0;
	Xgrid_test.xstep = xstep_test;
	Xgrid_test.fill_array();

	Array_1D <double> Pgrid_test(2000);
	Pgrid_test.xinit = 0.0;
	Pgrid_test.xstep = pstep_test;
	Pgrid_test.fill_array();

	for (ii = 2000/2; ii < 2000; ii++){
		Pgrid_test.grid[ii] = -1.0*Pgrid_test.grid[2000 - ii];
	}

	for (ii = 0; ii < 2000; ii++){
		Potential_1D_Test.grid[ii] = 0.5*Xgrid_test.grid[ii]*Xgrid_test.grid[ii];
		Kinetic_1D_Test.grid[ii] = H_Kinetic(Pgrid_test.grid[ii],0.0,Phys_Param);
	}
	EigenArray_1D <cplx> eigen_test(2000,10);
	eigen_test.xstep = xstep_test;

	get_wvfxns(eigen_test, Xgrid_test, Potential_1D_Test, Kinetic_1D_Test, Phys_Param, 1.0);

	ofstream ITPoutfile;
	ITPoutfile.open("output_data/ITP_Results.txt");
	for (kk = 0; kk < 2000; kk++)
	{
		ITPoutfile << Xgrid_test.grid[kk] << "\t" << real((eigen_test.get_array_addr(0))[kk]);
		ITPoutfile << "\t" << real((eigen_test.get_array_addr(1))[kk]) << "\t";
		ITPoutfile << real((eigen_test.get_array_addr(2))[kk]) << "\t";
		ITPoutfile << real((eigen_test.get_array_addr(3))[kk]) << "\t";
		ITPoutfile << real((eigen_test.get_array_addr(4))[kk]) << "\t";
		ITPoutfile << real((eigen_test.get_array_addr(5))[kk]) << "\t";
		ITPoutfile << real((eigen_test.get_array_addr(6))[kk]) << "\t";
		ITPoutfile << real((eigen_test.get_array_addr(7))[kk]) << "\t";
		ITPoutfile << real((eigen_test.get_array_addr(8))[kk]) << "\t";
		ITPoutfile << Pgrid_test.grid[kk] << "\t" << Potential_1D_Test.grid[kk];
		ITPoutfile << "\t" << Kinetic_1D_Test.grid[kk] << endl;
	}
	ITPoutfile.close();
*/

/*******************************************************************
	End Test Code Section
*********************************************************************/

	return 0;
}
