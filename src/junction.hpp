#ifndef SPLITOP_JUNCTION
#define SPLITOP_JUNCTION

#include "input_parser.hpp"
#include "wvfxn.hpp"

#include "junction.hpp"

using namespace std;

double HKinetic1D(double pp, double mass);
double HKinetic2D(double pp, double qq, double mass1, double mass2);
double VM (double ZZ, programInputs &IP);
double WL (double zz, programInputs &IP);
double WR (double zz, programInputs &IP);
double VVac (double zz, programInputs &IP);
double Ve (double zz, programInputs &IP);
double VeMS (double ZZ, programInputs &IP);
double Wcoup (double zz, double ZZ, programInputs &IP);
double ElecParab (double zz, double ZZ, programInputs &IP);
double VeM (double zz, double ZZ, programInputs &IP);
double V0 (double zz, double ZZ, programInputs &IP);
double absorbingPotential(double XX, programInputs &IP);
double absorbingPotentialMol(double XX, programInputs &IP);
double HPotentialNoAbs (double zz, double ZZ, programInputs &IP);
cplx HPotential (double zz, double ZZ, programInputs &IP);
double wvfxnNuclear (double ZZ, programInputs &IP);
cplx wvfxnElectron (double zz, programInputs &IP);
double hermite(int n, double x);
double factorial(int n);
double hoWvfxn(int n, double x, double xcen, double omega, double m, double hbar);
int kron_delta(int x, int y);


#endif