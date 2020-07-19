#include <sys/stat.h>
#include <cstdlib> 
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>
#include <mpi.h>
#include <Eigen/Dense>
#include "nsolver/nsolver.h"
#include "deDSI.h"

#include "cfbasics/cfvector.h"
#include <cfbasics/cfbasics.h>
using namespace std;

// This program calculates eigenvalues of fixed point 
// using Arnoldi iteration. The ideas and algorithm are based on Divakar
// Viswanath, "Recurrent motions within plane Couette turbulence",
// J. Fluid Mech.</em> 580 (2007), http://arxiv.org/abs/physics/0604062.

using namespace ded_ns;
using namespace nsolver;

int main( int argc, char *argv[] ){

    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm MPI_COMM_WORLD;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank != 0 ) { 
    std::ofstream* sink = new std::ofstream ("/dev/null");
    std::cout.rdbuf (sink->rdbuf());
    std::cerr.rdbuf (sink->rdbuf());
    } 

    string purpose("Compute spectrum of eigenvalues of equilibria, traveling waves, or periodic orbit using Arnoldi iteration");
    ArgList args(argc, argv, purpose);
    EigenvalsFlags eigenflags (args);

    Eigenvals E(eigenflags);

    unique_ptr<Newton> N;
    NewtonSearchFlags searchflags(args);
    searchflags.save(searchflags.outdir);
    N = unique_ptr<Newton>(new NewtonAlgorithm(searchflags));
    
    args.section("Program options");
    
    const string sigstr = args.getstr ("-sigma", "--sigma", "", "file containing sigma of sigma f^T(u) - u = 0 (default == identity)");
    
    const int  seed     = args.getint ("-sd",  "--seed",       1, "seed for random number generator");
    const Real EPS_du   = args.getreal("-edu","--epsdu", 1e-7, "magnitude of perturbation for numerical approximation of the Jacobian");
    const string duname = args.getstr ("-du", "--perturb", "", "initial perturbation field, random if unset");//perturbation of the jacobian
    
     const bool verbose         = args.getflag("-ver",      "--Verbose",  "make dsi verbose");
     const int Nx              = args.getint ("-Nx",     "--Xresolution", 512,  "Number of points to resolve the x-dimension");
     const int Ny              = args.getint ("-Ny",     "--Yresolution", 256,  "Number of points to resolve the y-dimension");
     const int D               = args.getint ("-Dim",    "--Dimension",     2,  "Dimension of the system (1d and 2d is implemented)");
     const Real Lx             = args.getreal("-valLx",     "--Xlength",       1,  "Length of the x-dimension");
     const Real Ly             = args.getreal("-valLy",     "--Ylength",       1,  "Length of the y-dimension");
     const Real Ra             = args.getreal("-valRa",     "--Rayleigh",  7.524e05, "Rayleigh number");
     const Real Be             = args.getreal("-valBe",     "--beta",  3e4, "beta");
     const Real Pr             = args.getreal("-valPr",     "--Prandtl",  1, "prandtl");
     const Real Cc             = args.getreal("-valCc",     "--C",  0, "a geometric parameter");
     const Real dtmax          = args.getreal("-Dtmax",  "--maxtimestep",  1e-4, "maximum time step");
     const Real dtmin          = args.getreal("-Dtmin",  "--minetimestep", 1e-7, "minimum time step");
     const Real timestep       = args.getreal("-tstep",  "--timestep"   , 1e-5, "initial time step");
     const Real T              = args.getreal("-T",      "--Time",  0.079, "Duration of a dns");
     const Real precondfPsi    = args.getreal("-pPsi",   "--PrecondfPsi",    1.0,  "preconditioning factor for Psi");
     const Real precondfTheta  = args.getreal("-pTheta", "--PrecondfTheta",  1.0,  "preconditioning factor for Theta");
     const Real precondfZeta   = args.getreal("-pZeta",  "--PrecondfZeta",   1.0,  "preconditioning factor for Zeta");
     const int solnum       = args.getint ("-Sn",     "--Solnum", -1, "solution number to read from file ");
     const string ifilename  = args.getstr ("-fname",  "--Filename",  "", " name and adress of the file to read the guess from"); 
     const string ofilepath = args.getstr ("-Opath",  "--OutputPath", "./", " output directory"); // TODO make it relative
     const string ofilename = args.getstr ("-Oname",  "--OutputFilename", "best", "output file base");
 
     const Real ax             = args.getreal("-ax",     "--xshift",  0.38432,   "initial guess for the shift in x ");
     const bool norm_max      = args.getflag ("-Nm",  "--Nmax" ,   "Normalize each field with respect to maximum of that field");
     const bool norm_l2       = args.getflag ("-Nl",  "--Nl2norm" ,"Normalize each field with respect to the l2norm of that field");
     const Real amp           = args.getreal("-A", "--amp", 1, "factor of rescaling");
    
    args.check();
    args.save("./");
    args.save (eigenflags.outdir);

     
    bool Rxsearch = searchflags.xrelative;
    bool Tsearch  = searchflags.solntype == PeriodicOrbit ? true : false;
  
    if(Rxsearch && Tsearch == true) cout << "newton: ax is " << ax << " T (orbital) is " << T << endl; 
    cout << "Rxsearch is " << Rxsearch << " Tsearch is " << Tsearch << endl; 
     
    DEDSI dsi(Rxsearch, Tsearch, norm_max, norm_l2, ax);
    dsi.verbose(verbose);
   
    dsi.setSystem( Ra, Be, Pr, Cc , Nx, Ny, Lx, Ly, D, T , dtmax, dtmin, timestep, precondfPsi, precondfTheta, precondfZeta, ofilepath, ofilename);
    bp::object soln = dsi.initialGuess(ifilename, solnum, amp);
    VectorXd xSoln = dsi.makeVector(soln);

    srand48(seed);


    // Set up DNS operator ("A" in Arnoldi A*b terms)
    if(rank == 0)
      cout << "setting up DNS and initial fields..." << endl;
    
    const Real l2u = dsi.DSIL2Norm(xSoln); //FIXME #1 how do I call L2Norm?

    const Real eps = ( l2u < EPS_du ) ? EPS_du: EPS_du/l2u; //for choice of epsilon, see eq. (15) in reference
					//C.J. Mack, P.J. Schmid/Journal of Computational Physics 229 (2010) 541-560
    if(rank == 0)
      cout << "computing sigma f^T(u)..." << endl;
    
    //Construct the dynamical-systems interface object depending on the given parameters.
    
    // Check if sigma f^T(u) - u = 0
    
    VectorXd Gx = dsi.eval(xSoln);
    
    Real l2normGx = chflow::L2Norm(Gx);

    if(rank == 0) 
      cout << "L2Norm(u - sigma f^T(u)) = " << l2normGx  << endl;
    
    if (l2normGx > 1e-06)
     cferror("error: (u, sigma, T) is not a solution such as sigma f^T(u) - u = 0");
    
    
    VectorXd dx = VectorXd::Random(xSoln.size()); 

    printout("L2Norm(du) = " + r2s(dsi.DSIL2Norm(dx)));
    printout("rescaling du by eps_du = " + r2s(EPS_du));
    dx *= EPS_du/dsi.DSIL2Norm(dx); //FIXME rescaling must be by integration l2norm so has to be physical
    printout("L2Norm(du) = " + r2s(dsi.DSIL2Norm(dx)));
  
    
    E.solve(dsi, xSoln, dx, T, eps);
    
    MPI_Finalize();

return 0;

}

