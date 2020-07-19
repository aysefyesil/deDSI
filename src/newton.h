#pragma once
#include "nsolver/nsolver.h"
#include "deDSI.h"

namespace ded_ns {

/** \brief Perform Newton search with nsolver Newton-Krylov-Hookstep methdo
 * \param[in,out] w initial guess for fixed point
 * \param[in] C Cylinder describing the problem and equations
 * \param[in] hflags search flags (see nsolver doc)
 * \param[in] iflags Integration flags, used for computing w(T) - w(0) as the function to be solved
 * \param[in] T as in w(T) - w(0). If T == 0, use C.rhs(w) instead.
 */
// here integer d is the dimension of the system coming from the Dedalus program	
//
Real newtonHookstepSearch (chflow::NewtonSearchFlags searchflags, double T, double dtmax, double dtmin, double timestep, string log, int Nx, int Ny, Real Lx, Real Ly, int d, Real Ra, Real Be, Real Pr, Real Cc, double precondfPsi, double precondfTheta, double precondfZeta, int solnum , string ifilename, string ofilepath, string ofilename, int rank, bool norm_max, bool norm_l2, double ax, bool verbose, Real amp) ;

Real continuation (chflow::NewtonSearchFlags searchflags, ContinuationFlags cflags, double T, double dtmax, double dtmin, double timestep, int Nx, int Ny, string muname, Real Lx, Real Ly, int d, Real Ra, Real Be, Real Pr, Real Cc, double precondfPsi, double precondfTheta, double precondfZeta, int solnum , string ifilename, string ofilepath, string ofilename, int rank, bool norm_max, bool norm_l2, double ax, bool Forward, bool Backward, bool verbose, Real amp);



} // namespace
