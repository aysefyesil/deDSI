#pragma once
#include "nsolver/nsolver.h"
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <Python.h>
#include "error.h"
#include <fstream>

namespace bp = boost::python;
namespace np = boost::python::numpy;

namespace ded_ns {


enum class continuationParameter
{ LX, LY, RA, BE, PR, CC };

  
//class DEDSI : public nsolver::DSI {
class DEDSI : public DSI {
public:

  DEDSI();

  DEDSI(bool xrelative, bool Tsearch, bool norm_max, bool norm_l2, double ax);

  virtual ~DEDSI(){}; 
 

  ofstream fout1;
  ofstream fout2;
  ofstream fout3;

  void display(np::ndarray &u);
  void display(bp::object &u);
  void verbose(bool v);

  void setSystem(Real Ra, Real Be, Real Pr, Real Cc, int Nx, int Ny, Real Lx, Real Ly, int d, double T , double dtmax, double dtmin, double timestep, double precondfPsi, double precondfTheta, double precondfZeta, string outpath, string outfile);

  bp::object initialGuess(const string filename, int read_it_from_file, const Real Amp );

  VectorXd precondition(const VectorXd& x);
  VectorXd unPrecondition(const VectorXd& x);

  np::ndarray vectorXdToNumpyArray1d( const VectorXd & v );
  np::ndarray vectorXdToNumpyArray2d( const VectorXd & v );
  np::ndarray vectorXdToNumpyArray( const VectorXd & v );

  VectorXd numpyArrayToVectorXd1d( const np::ndarray & na );
  VectorXd numpyArrayToVectorXd2d( const np::ndarray & na );
  VectorXd numpyArrayToVectorXd( const np::ndarray & na );

  VectorXd boostObjectToVectorXd1d( const bp::object & bpo );
  VectorXd boostObjectToVectorXd2d( const bp::object & bpo );
  VectorXd boostObjectToVectorXd( const bp::object & bpo );

  
  void buildSolver( int Nx, int Ny, Real Lx, Real Ly , Real Ra, Real Be, Real Pr, Real Cc);
  VectorXd normalize( VectorXd & x, bool norm_max, bool norm_l2, double pnorm, double tnorm, double znorm);
  VectorXd unNormalize( np::ndarray & u, double pnorm, double tnorm, double znorm );

  bp::object translator(const np::ndarray & a); 
  bp::object reflector(const np::ndarray & a); 

  void setLocalInfo();
  int rank();

  int getSize();

  void setTime( double T, double dtmax, double dtmin, double timestep );
  void setSigma( int sp, int st, int sx, int sy, double ax, double Torb );
  void setPrecondf( double precondfPsi, double precondfTheta, double precondfZeta);
  
  void adjust_for_T (double T);

  VectorXd makeVector ( const bp::object & u );
  void extractVector (const VectorXd& x); 
  np::ndarray makeField( const VectorXd& x );

  VectorXd eval ( const VectorXd& x );
  VectorXd G ( const VectorXd& x );
  VectorXd f ( const np::ndarray& u ,  double T);


  VectorXd sigmaOps ( const VectorXd& fu );

//void save ( const VectorXd& x, const string filebase, const string outdir ) ;
  void save (const VectorXd & x, const std::string filebase, const std::string outdir, const bool fieldsonly);
  void saveEigenvec (const VectorXd& x, const string label, const string outdir) ;
  void saveEigenvec (const VectorXd& x1, const VectorXd& x2, const string label1, const string label2, const string outdir) ;

  void initiateSave(string foldername ,const string filename);

  bp::object otp; //pointer to the file output

  Real observable( VectorXd& x );
  Real observable2( VectorXd& x );
  Real DSIL2Norm( const VectorXd& x);
  Real DSIL2Norm_X( const VectorXd& x); //TODO decide if this stays or not

  string statsHeader() ;
  string stats( const VectorXd& x );

  VectorXd xdiff ( const VectorXd& u ) ;
  VectorXd tdiff ( const VectorXd& a , Real epsDt ) ;


  Real getMu();
  Real mu() const ;
  void setMu( Real mu );
  void updateMu( Real mu );

  void chooseMu (continuationParameter mu) ;
  continuationParameter s2cPar (string muname) ;
  void chooseMu (string muName) ;


  bp::object deModule ; 


  bp::object setSolver;
  bp::object solver ;

  bp::object integrate;
  bp::object uT;

  bp::object exact;
  bp::object xT_exact; 
  
  bp::object setparamMu;

  bp::object xDiff;
  bp::object dud_x; 

  bp::object translate;
  bp::object reflect;
  bp::object sfu ;

  int rank_ =0;

  double pnorm_ = 1;
  double tnorm_ = 1;
  double znorm_ = 1;

protected:
  bool normal_ = true;
  bool norm_max_ ;
  bool norm_l2_ ;

  int knt_ = 1;

  double T_ ;
  double dt_;
  double dtmax_;
  double dtmin_;

  int Nx_ = 1; 
  int Ny_ = 1;
  int Ny_p_ = 1;
  

  int kntr_ = 0;
  int kntr_vxdtna_ = 0;

  Real Lx_;
  Real Ly_;
  Real Ra_; 
  Real Be_; 
  Real Pr_; 
  Real Cc_; 


  int d_ ;

  int n_ = 0;
  int ss_ = 0;
  int two_fields_ = 1;
  int three_fields_ = 1;

  int sp_ = 1;
  int st_ = 1;
  int sx_ = 1;
  int sy_ = 1;
 
  //bool precondf_ = true;  
  bool precondf_ = false;  
  double precondfPsi_   = 1;
  double precondfTheta_ = 1;
  double precondfZeta_  = 1;


  bool xrelative_ ;
  double ax_ = 0.;

  bool Tsearch_;
  bool verbose_ = false;
  bool normalize_ = false;

  VectorXd silX_ ;
  int kauntir = 0;
  continuationParameter cPar_ = continuationParameter::RA;
};

}// namespace ded_ns 

