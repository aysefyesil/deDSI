#include "newton.h"
#include "deDSI.h"
#include <cfbasics/cfbasics.h>

using namespace nsolver;

namespace ded_ns {

Real newtonHookstepSearch (chflow::NewtonSearchFlags searchflags, double T, double dtmax, double dtmin, double timestep, string log, int Nx, int Ny, Real Lx, Real Ly, int d, Real Ra, Real Be, Real Pr, Real Cc, double precondfPsi, double precondfTheta, double precondfZeta, int solnum , string ifilename, string ofilepath, string ofilename, int rank, bool norm_max, bool norm_l2, double ax, bool verbose, Real amp) {
  ostream* out;
  if (log == "none") {
    out = new ofstream("/dev/null");
  } else if (log == "cout")
    out = &cout;
  else
    out = new ofstream(log);
  searchflags.logstream = out;

  bool Rxsearch = searchflags.xrelative;
  bool Tsearch  = searchflags.solntype == PeriodicOrbit ? true : false;

  if(Rxsearch && Tsearch == true) cout << "newton: ax is " << ax << " T (orbital) is " << T << endl; 
  cout << "Rxsearch is " << Rxsearch << " Tsearch is " << Tsearch << endl; 
   
  DEDSI dsi(Rxsearch, Tsearch, norm_max, norm_l2, ax);
  dsi.verbose(verbose);
 
  dsi.setSystem( Ra, Be, Pr, Cc , Nx, Ny, Lx, Ly, d, T , dtmax, dtmin, timestep, precondfPsi, precondfTheta, precondfZeta, ofilepath, ofilename);
  bp::object guess = dsi.initialGuess(ifilename, solnum, amp);

  VectorXd xGuess = dsi.makeVector(guess);

  cout << std::scientific << std::setprecision(17) << "DSIL2Norm of the guess l2norm(x) " << chflow::L2Norm(xGuess) << " l2norm(u) " << dsi.DSIL2Norm(xGuess) << " x(0) " << xGuess(0) <<endl; 


  unique_ptr<Newton> N;
  N = unique_ptr<Newton>(new NewtonAlgorithm(searchflags));

  Real residual = 0;

  cout << "in newton: solve is started " << Tsearch << endl; 

  N->solve(dsi, xGuess, residual);

  return residual;
}


Real continuation (chflow::NewtonSearchFlags searchflags, ContinuationFlags cflags, double T, double dtmax, double dtmin, double timestep, int Nx, int Ny, string muname, Real Lx, Real Ly, int d, Real Ra, Real Be, Real Pr, Real Cc, double precondfPsi, double precondfTheta, double precondfZeta, int solnum , string ifilename, string ofilepath, string ofilename, int rank, bool norm_max, bool norm_l2, double ax, bool Forward, bool Backward, bool verbose, Real amp){

  if(rank==0) cout << "Starting continuation" << endl;

  cfarray<Real> mu(3);
  cfarray<VectorXd> x(3);


  bool Rxsearch = searchflags.xrelative;
  bool Tsearch  = searchflags.solntype == PeriodicOrbit ? true : false;

  cout << "continuation: ax is " << ax << " T (orbital) is " << T << " Mu (Rayleigh) is " << Ra << endl; 

  DEDSI dsi(Rxsearch, Tsearch, norm_max, norm_l2, ax);

  dsi.setSystem( Ra, Be, Pr, Cc , Nx, Ny, Lx, Ly, d, T , dtmax, dtmin, timestep, precondfPsi, precondfTheta, precondfZeta, ofilepath, ofilename);

  bp::object guess = dsi.initialGuess(ifilename, solnum, amp);
  

  VectorXd xGuess = dsi.makeVector(guess);
  cout << " continuation size of xGuess is " << xGuess.size() << endl;

  unique_ptr<Newton> N;
  N = unique_ptr<Newton>(new NewtonAlgorithm(searchflags));


  Real dmu = cflags.initialParamStep;
  dsi.chooseMu(muname);
  Real mu0 = dsi.mu();
  cout << " newton.cpp mu0: " << mu0 << endl;
  
 
  VectorXd x0 = xGuess;

  for (int i=0; i<3; ++i) {

     int step = i - 1;

     mu[i] = mu0 - step*dmu;

    if(Forward && i == 2){//ascending
        mu[0] = mu0 + 2 * dmu ; 
        mu[1] = mu0 + dmu ;
        mu[2] = mu0; 
    }
    else if(Backward && i == 2){//descending
        mu[0] = mu0  - 2 * dmu ; 
        mu[1] = mu0 - dmu ;
        mu[2] = mu0 ;
    }

    x[i] = x0;
  }
  
 return continuation(dsi, *N, x, mu, cflags);

} 


} // namespace 

