typedef std::complex<double> cplx;

//Kinetic energy operator for the total Hamiltonian in momentum space
double H_Kinetic(double pp, double qq, Phys_Parameters &P);

//Potential energy function for the molecule, Lennard-Jones 6-12 potential
double V_M (double ZZ, Phys_Parameters &P);

//Junction turn-on function
double W_L (double zz, Phys_Parameters &P);

//Junction turn-off function
double W_R (double zz, Phys_Parameters &P);

//Vacuum potential function
double V_Vac (double zz, Phys_Parameters &P);

//Electronic barrier potential energy function
double V_e (double zz, Phys_Parameters &P);

//Electron-Molecule-Surface interaction via the Antoniewicz mechanism
double V_eMS (double ZZ, Phys_Parameters &P);

//Coupling function, either Gaussian or Parabolic, dependent on P.coupling parameter
double W_coup (double zz, double ZZ, Phys_Parameters &P);

//Parabolic coupling function without turn-on and turn-off functions
double Elec_parab (double zz, double ZZ, Phys_Parameters &P);

//Electron-Molecule Coupling term
double V_eM (double zz, double ZZ, Phys_Parameters &P);

//Electron and Molecular Potential terms
double V_0 (double zz, double ZZ, Phys_Parameters &P);

//Imaginary absorbing potentials for electron, molecule, Manolopoulos2002 and Manolopoulos2004
double absorbing_potential(double XX, Phys_Parameters &P, Comp_Parameters &C);

double absorbing_potential_mol(double XX, Comp_Parameters &C);

//Full Hamiltonian potential term, without the absorbing potential
double H_Potential_NoAbs (double zz, double ZZ, Phys_Parameters &P);

//Full Hamiltonian potential term, with absorbing potential
cplx H_Potential (double zz, double ZZ, Phys_Parameters &P, Comp_Parameters &C);

//Returns ground state harmonic oscillator function for the phonon mode
double wvfxn_nuclear (double ZZ, Phys_Parameters &P);

//Returns an electron wavepacket
cplx wvfxn_electron (double zz, Phys_Parameters &P);
