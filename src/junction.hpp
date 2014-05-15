/**
 * \file junction.hpp
 * \brief Junction specific routines
 *
 * \author Joshua E.\ Szekely
 * \date May 2014
 */

#ifndef SPLITOP_JUNCTION
#define SPLITOP_JUNCTION

#include "input_parser.hpp"
#include "wvfxn.hpp"

#include "junction.hpp"

using namespace std;

/**
 * @brief 1D kinetic energy operator
 * @details Basic kinetic energy operator, \f$ \frac{\hbar^2}{2m} \frac{\partial^2}{\partial z^2}\f$
 *
 * @param pp momentum
 * @param mass particle mass
 *
 * @return kinetic energy
 */
double HKinetic1D(double pp, double mass);

/**
 * @brief 2D kinetic energy operator
 * @details Two dimensional kinetic energy operator, \f$ \frac{\hbar^2}{2m_1} \frac{\partial^2}{\partial z^2} + \frac{\hbar^2}{2m_2} \frac{\partial^2}{\partial Z^2}\f$
 *
 * @param pp momentum of dimension 1
 * @param qq momentum of dimension 2
 * @param mass1 mass in dimension 1
 * @param mass2 mass in dimension 2
 * @return kinetic energy
 */
double HKinetic2D(double pp, double qq, double mass1, double mass2);

/**
 * @brief Molecular potential energy surface
 * @details The molecular/phonon potential energy surface described by a
 * Lennard-Jones 6-12 potential
 *
 * @param ZZ Position
 * @param IP Input parameters describing the well shape
 *
 * @return Potential energy
 */
double VM (double ZZ, programInputs &IP);
/**
 * @brief Smooth turn on function
 * @details \f$w_L(z) = (1+e^{\alpha(z-z_L)})^{-1}\f$
 *
 * @param zz Position
 * @param IP Input parameters
 *
 * @return Potential energy
 */
double WL (double zz, programInputs &IP);
/**
 * @brief Smooth turn off function
 * @details \f$w_L(z) = (1+e^{\alpha(z-z_L)})^{-1}\f$
 *
 * @param zz Position
 * @param IP Input Parameters
 *
 * @return Potential Energy
 */
double WR (double zz, programInputs &IP);
/**
 * @brief Vacuum potential energy
 * @details \f$V_{Vac} = eV_b(1-(z-z_L)/L)\f$
 *
 * @param zz Position
 * @param IP Input Parameters
 *
 * @return Potential Energy
 */
double VVac (double zz, programInputs &IP);
/**
 * @brief Electron potential energy
 * @details \f$V_e(z) = w_R(z)V_R + w_L(z)V_L + [1-w_R(z)-w_L(z)]V_{Vac}(z)\f$
 *
 * @param zz Position
 * @param IP Input Parameters
 *
 * @return Potential Energy
 */
double Ve (double zz, programInputs &IP);
/**
 * @brief Electron Molecule Surface Coupling
 * @details \f$ V_{e-p-s}(Z) = E_A - \frac{C_{coup}}{Z-z_L} \f$
 *
 * @param ZZ Position
 * @param IP Input Parameters
 *
 * @return Coupling Potential Energyu
 */
double VeMS (double ZZ, programInputs &IP);
/**
 * @brief Coupling term describing the potential energy well shape
 * @details \f{eqnarray*}{
 * 			W(z,Z) 	&=& e^{\frac{(z-Z)^2}{2\sigma^2}} \\
 * 					&=& w_{on}(z) w_{off}(z) \left(1-\left(\frac{z-Z}{r_{rad}}\right)^2\right)
 *			\f}
 *
 * @param zz Electronic Position
 * @param ZZ Molecular Position
 * @param IP Input parameters
 * @return Weighting Factor
 */
double Wcoup (double zz, double ZZ, programInputs &IP);
/**
 * @brief Parabolic electron potential
 * @details Quadtratic potential well for 1D calculation of electronic eigenstates
 * Used with imaginary time propagation routines
 *
 * @param zz Electron Position
 * @param ZZ Molecule Position
 * @param IP Input Parameters
 * @return Potential Energy
 */
double ElecParab (double zz, double ZZ, programInputs &IP);
/**
 * @brief Electron Molecule coupling term
 * @details  \f$V_{e-p}(z,Z) = W(z,Z)\left[\hat V_{e-p-s}(Z) - \hat V_e(z) \right]\f$
 *
 * @param zz Electron Position
 * @param ZZ Molecule Position
 * @param IP Input Parameters
 * @return Potential Energy
 */
double VeM (double zz, double ZZ, programInputs &IP);
/**
 * @brief Potential surface without coupling
 * @details \f$ V_0(z,Z) = V_e(z) + V_M(Z)\f$
 *
 * @param zz Electron Position
 * @param ZZ Molecule Position
 * @param IP Input Parameters
 * @return Potential Energy
 */
double V0 (double zz, double ZZ, programInputs &IP);
/**
 * @brief Absorbing boundary conditions tuned for the electronic subspace
 * @details See ManolopoulosJCP2002 for detatils
 *
 * @param XX Position
 * @param IP Input Parameters
 *
 * @return Absorbing Potential
 */
double absorbingPotential(double XX, programInputs &IP);
/**
 * @brief Absorbing boundary conditions tuned for the molecular subspace
 * @details See ManolopoulosJCP2002 for detatils
 *
 * @param XX Position
 * @param IP Input Parameters
 *
 * @return Absorbing Potential
 */
double absorbingPotentialMol(double XX, programInputs &IP);
/**
 * @brief Full Potential Energy Junction Hamiltonian, no absorbing boundary
 *
 * @param zz Electron Position
 * @param ZZ Molecule Position
 * @param IP Input Parameters
 * @return Potential Energy Surface
 */
double HPotentialNoAbs (double zz, double ZZ, programInputs &IP);
/**
 * @brief Full Potential Energy Junction Hamiltonian, with absorbing boundary
 *
 * @param zz Electron Position
 * @param ZZ Molecule Position
 * @param IP Input Parameters
 * @return Potential Energy Surface
 */
cplx   HPotential (double zz, double ZZ, programInputs &IP);
/**
 * @brief Calculates an approximate molecular wavefunction
 *
 * @param ZZ Molecule Position
 * @param IP Input Parameters
 *
 * @return complex valued wavefunction
 */
double wvfxnNuclear (double ZZ, programInputs &IP);
/**
 * @brief Electron wavepacket
 * @details \f$ \psi(z,Z;t=0) = N e^{-\frac{i}{\hbar}p_0 z - \frac{(z-z_0)^2}{2\xi^2}} \phi_{molec}^0(Z) \f$
 *
 * @param zz Electron Position
 * @param IP Input Parameters
 *
 * @return complex valued wavefunction
 */
cplx   wvfxnElectron (double zz, programInputs &IP);
/**
 * @brief Hermite Polynomials
 * @details Recursively generated Hermite polynomials for harmonic oscillator wavefunctions
 *
 * @param n Polynomial degree
 * @param x value
 *
 * @return value of the nth Hermite polynomial evaluated at x
 */
double hermite(int n, double x);
/**
 * @brief Factorial function
 *
 * @return n!
 */
double factorial(int n);
/**
 * @brief Harmonic oscillator wavefunctions
 *
 * @param n target eigenstate
 * @param x position
 * @param xcen center of harmonic potential
 * @param omega oscillator frequency
 * @param m mass
 * @param hbar \f$ \hbar \f$
 * @return nth harmonic oscillator eigenfunction evaluated at x
 */
double hoWvfxn(int n, double x, double xcen, double omega, double m, double hbar);
/**
 * @brief Kronecker Delta function
 * @details \f{eqnarray*}{
 * 			\delta_{xy}	&=& 1 if x = y \\
 * 						&=& 0 if x \ne y
 * 			\f}
 *
 * @param x first integer
 * @param y second integer
 *
 * @return 0 or 1
 */
int    kron_delta(int x, int y);


#endif