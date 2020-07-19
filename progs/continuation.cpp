#include "../src/newton.h"
//#include "../src/test.h"
#include <cstdlib> 
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>
#include <mpi.h>

#include <cfbasics/cfbasics.h>

using namespace ded_ns;

int main( int argc, char *argv[] ){

 MPI_Init(&argc, &argv);
 int rank;
 MPI_Comm MPI_COMM_WORLD;
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank != 0 ) { 
    // Open ofstream on heap to preserve it after return
    std::ofstream* sink = new std::ofstream ("/dev/null");

    // Mute standard output
    std::cout.rdbuf (sink->rdbuf());
    // Optionally mute standard error
    std::cerr.rdbuf (sink->rdbuf());
  } 
 
    ArgList args (argc, argv, "Parameters for Busse Annulus");


    ContinuationFlags cflags(args);
    cflags.save();

    NewtonSearchFlags searchflags(args);
    searchflags.save(searchflags.outdir);
 

     const bool verbose         = args.getflag("-ver",      "--Verbose",  "make dsi verbose");
  
     const bool Forward        = args.getflag("-Fwd",    "--Forward",       "Continue the parameter forward");
     const bool Backward       = args.getflag("-Bwd",    "--Backward",      "Continue the parameter backward");

     const int Nx              = args.getint ("-Nx",     "--Xresolution", 512,  "Number of points to resolve the x-dimension");
     const int Ny              = args.getint ("-Ny",     "--Yresolution", 256,  "Number of points to resolve the y-dimension");
     const int D               = args.getint ("-Dim",    "--Dimension",     2,  "Dimension of the system (1d and 2d is implemented)");
     
     const string muname       = args.getstr("-cont","--continuation","", "Parameters can be: LX(size of dim. x), LY(size of dim. y), RA (Rayleigh), BE (beta), PR (Prandtl), CC (C)");
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

     const Real ax             = args.getreal("-ax",     "--xshift",  0.38432,   "initial guess for the shift in x ");

     int solnum       = args.getint ("-Sn",     "--Solnum", -1, "solution number to read from file ");
     string filename  = args.getstr ("-fname",  "--Filename",  "", " name and adress of the file to read the guess from"); 
     string ofilepath = args.getstr ("-Opath",  "--OutputPath", "./", " output directory"); 
     string ofilename = args.getstr ("-Oname",  "--OutputFilename", "best", "output file base");

     const bool norm_max      = args.getflag ("-Nm",  "--Nmax" ,   "Normalize each field with respect to maximum of that field");
     const bool norm_l2       = args.getflag ("-Nl",  "--Nl2norm" ,"Normalize each field with respect to the l2norm of that field");
     const Real amp           = args.getreal("-A", "--amp", 1, "factor of rescaling");

 args.check();
 args.save();

 cout << "Starting point: " << endl;
 cout << "Nx: " << Nx << " Ny: " << Ny << " Ra: " << Ra << " beta: " << Be << " Prandtl: " << Pr << " C: " << Cc << " precondtioning factors for psi: " << precondfPsi << " for theta: " << precondfTheta << " for zeta: " << precondfZeta << endl; 
 cout << "T : " << T << endl;
 cout << "Reading from the file: " << filename << " the solution number for guess is " << solnum << endl;

 string log = "cout";
 
 continuation ( searchflags, cflags, T, dtmax, dtmin, timestep, Nx, Ny, muname, Lx, Ly, D, Ra, Be, Pr, Cc, precondfPsi, precondfTheta, precondfZeta, solnum , filename, ofilepath, ofilename, rank, norm_max, norm_l2, ax, Forward, Backward, verbose, amp);
 
 MPI_Finalize();

return 0;
}
