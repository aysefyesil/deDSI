#include "deDSI.h"
#include <math.h>
using namespace std;


namespace ded_ns {
// the functions that are directly communicating with the nsolver should always return 
// vector field  

DEDSI::DEDSI()
{

  try{
     Py_Initialize();
     np::initialize();
     bp::import("sys").attr("path").attr("append")("/home/yesil/deDSI_SRC/ded/"); // add location of de_module_x to pythonpath, otherwise compiler searches the file in the python path
     } catch(bp::error_already_set &) {
          handle_pyerror();
     }
  try{
      deModule = bp::import("de_module_ns_dsi");
     } catch(bp::error_already_set &) {
          handle_pyerror();
     }
}

DEDSI::DEDSI(bool xrelative, bool Tsearch,bool norm_max, bool norm_l2, double ax)
 :
  norm_max_(norm_max),
  norm_l2_(norm_l2),
  xrelative_(xrelative),
  ax_(ax),
  Tsearch_(Tsearch)
{
  if(Tsearch_)   fout1.open("Tbest.asc");
  if(xrelative_) fout2.open("ax_best.asc");

  try{
     Py_Initialize();
     np::initialize();
     bp::import("sys").attr("path").attr("append")("/home/yesil/deDSI_SRC/ded/"); // add location of de_module_x to pythonpath, otherwise compiler searches the file in the python path
     } catch(bp::error_already_set &) {
          handle_pyerror();
     }
  try{
      deModule = bp::import("de_module_ns_dsi");
     } catch(bp::error_already_set &) {
          handle_pyerror();
     }
}


void DEDSI::verbose(bool v){
  verbose_ = v;
}

void DEDSI::setSystem(Real Ra, Real Be, Real Pr, Real Cc, int Nx, int Ny, Real Lx, Real Ly, int d, double T , double dtmax, double dtmin, double timestep, double precondfPsi, double precondfTheta, double precondfZeta, string outpath, string outfile){

  if(verbose_ ==false) cout << " This is set system " << endl; 
  d_ = d;                                                //  set the dimension in the beginning as well
  buildSolver(Nx, Ny, Lx, Ly, Ra, Be, Pr, Cc);
  
  setLocalInfo();                                        //this is necessary to set the local info
  setTime( T , dtmax, dtmin, timestep);                                   //before integrating this information should be set
  setPrecondf(precondfPsi, precondfTheta, precondfZeta); // preconditioning has to set
}


void DEDSI::buildSolver( int Nx, int Ny, Real Lx, Real Ly, Real Ra, Real Be, Real Pr, Real Cc){

  if(verbose_ == true) cout << " This is build solver " << endl; 
  
  setSolver = deModule.attr("build_solver");
  solver = setSolver(Nx , Ny, Lx, Ly, Ra, Be, Pr, Cc);

  Nx_ = Nx;
  Ny_ = Ny;

  Lx_ = Lx;
  Ly_ = Ly;
  Ra_ = Ra;
  Be_ = Be;
  Pr_ = Pr;
  Cc_ = Cc;

  Ny_p_ = Ny; //this is Ny per core. It will be reset when setSystem calls setLocalInfo() 
 
}


void DEDSI::initiateSave(string foldername, const string filename){

  if(verbose_ == true) cout << " This is initiate save " << endl; 
  bp::object in_wr_guess = deModule.attr("initiate_write_guess");
  otp = in_wr_guess(foldername, filename );
}

void DEDSI::setLocalInfo(){

//this function sets the local information about rank and Ny per core. 


  bp::object getlocinfo = deModule.attr("get_local_info");
  bp::object gli = getlocinfo(solver);

  rank_ = bp::extract<int>(gli[0]);
  if(verbose_ == true) cout << " This is set local info " << endl; 
  Ny_p_ = bp::extract<int>(gli[2]);

}	

int DEDSI::rank(){
 return rank_;
}

bp::object DEDSI::initialGuess(const string filename, int read_it_from_file, const Real Amp){


  if(verbose_ == true) cout << " This is set initial guess " << endl; 

  bp::object distribute_guess = deModule.attr("distribute_guess");

  bp::object guess = distribute_guess(solver, filename, read_it_from_file , Amp);

  return guess;
}

VectorXd DEDSI::normalize( VectorXd & x, bool norm_max, bool norm_l2, double pnorm, double tnorm, double znorm ){
     
  if(verbose_ == true) cout << " deDSI normalize " << endl;

  np::ndarray u = makeField(x);

  bp::object normalise = deModule.attr("normalization");
  bp::object normalized = normalise(solver, u , norm_max, norm_l2, pnorm, tnorm, znorm );
  VectorXd normal_x = boostObjectToVectorXd(normalized[0]);
  
  pnorm_ = bp::extract<double>(normalized[1]);
  tnorm_ = bp::extract<double>(normalized[2]);
  znorm_ = bp::extract<double>(normalized[3]);

  return normal_x;  
}

VectorXd DEDSI::unNormalize( np::ndarray & u, double pnorm, double tnorm, double znorm ){
     
  if(verbose_ == true ) cout << " deDSI unNormalize " << endl;

  bp::object unnormalize = deModule.attr("unnormalization");
  bp::object unnormalized = unnormalize(solver, u , pnorm, tnorm, znorm );
  VectorXd unnormal_x = boostObjectToVectorXd(unnormalized);

  return unnormal_x;  

}

Real DEDSI::mu() const{

  Real m ;

  if (cPar_ == continuationParameter::LX) {
     m = Lx_;
  } else if (cPar_ == continuationParameter::LY) {
     m = Ly_;
  } else if (cPar_ == continuationParameter::RA) {
     m = Ra_;  
  } else if (cPar_ == continuationParameter::BE) {
     m = Be_;
  } else if (cPar_ == continuationParameter::PR) {
     m = Pr_;
  } else if (cPar_ == continuationParameter::CC) {
     m = Cc_;
  } else {
    throw invalid_argument ("deDSI::updateMu(): continuation parameter is unknown");
  }

  if(verbose_ == true) cout << " This is mu() " << endl; 

  return m; 
}

void DEDSI::setMu(Real mu){

  if(verbose_ == true) cout << " This is set mu() " << endl; 

  setparamMu = deModule.attr("set_mu");
  setparamMu(mu);

  if (cPar_ == continuationParameter::LX) {
     Lx_ = mu;
  } else if (cPar_ == continuationParameter::LY) {
    Ly_ = mu;
  } else if (cPar_ == continuationParameter::RA) {
    Ra_ = mu;
  } else if (cPar_ == continuationParameter::BE) {
    Be_ = mu;
  } else if (cPar_ == continuationParameter::PR) {
    Pr_  =  mu;
  } else if (cPar_ == continuationParameter::CC) {
    Cc_ = mu;
  } else {
    throw invalid_argument ("deDSI::updateMu(): continuation parameter is unknown");
  }

}

int DEDSI::getSize(){
 
 int S = Nx_ * Ny_p_;
 if(two_fields_)   S = 2 * Nx_ * Ny_p_;
 if(three_fields_) S = 3 * Nx_ * Ny_p_;
 if(verbose_ == true) 
    cout << " deDSI: getSize() is called. Size of any field is " << S << endl;
 return S;
}

void DEDSI::updateMu (Real mu) {
  if(verbose_ == true ) 
     cout << "deDSI: updateMu is called" << endl;

  if (cPar_ == continuationParameter::LX) {
     Lx_ = mu;
  } else if (cPar_ == continuationParameter::LY) {
    Ly_ = mu;
  } else if (cPar_ == continuationParameter::RA) {
    Ra_ = mu;
  } else if (cPar_ == continuationParameter::BE) {
    Be_ = mu;
  } else if (cPar_ == continuationParameter::PR) {
    Pr_  =  mu;
  } else if (cPar_ == continuationParameter::CC) {
    Cc_ = mu;
  } else {
    throw invalid_argument ("deDSI::updateMu(): continuation parameter is unknown");
  }

  solver = setSolver(Nx_, Ny_, Lx_, Ly_, Ra_, Be_, Pr_, Cc_);
  if(verbose_ == true )
     cout <<"deDSI: Mu is " << Ra_ << endl;
}


void DEDSI::chooseMu (string muName) {
  chooseMu (s2cPar (muName));
}

void DEDSI::chooseMu (continuationParameter mu) {
  cPar_ = mu;
  switch (mu) {
    case continuationParameter::LX:
      updateMu (Lx_);
      break;
    case continuationParameter::LY:
      updateMu (Ly_);
      break;
    case continuationParameter::RA:
      updateMu (Ra_);
      break;
    case continuationParameter::BE:
      updateMu (Be_);
      break;
    case continuationParameter::PR:
      updateMu (Pr_);
      break;
    case continuationParameter::CC:
      updateMu (Cc_);
      break;
    default:
      throw invalid_argument ("deDSI::chooseMu(): continuation parameter is unknown");
  }
}


continuationParameter DEDSI::s2cPar (string muname) {
  std::transform (muname.begin(), muname.end(), muname.begin(), ::tolower);
  if (muname == "lx")
    return continuationParameter::LX;
  else if (muname == "ly")
    return continuationParameter::LY;
  else if (muname == "ra")
    return continuationParameter::RA;
  else if (muname == "be")
    return continuationParameter::BE;
  else if (muname == "pr")
    return continuationParameter::PR;
  else if (muname == "cc")
    return continuationParameter::CC;
  else
    throw invalid_argument ("deDSI::s2cPar(): continuation parameter '"+muname+"' is unknown");
}



void DEDSI::setPrecondf( double precondfPsi, double precondfTheta, double precondfZeta){
  precondfPsi_   = precondfPsi;
  precondfTheta_ = precondfTheta;
  precondfZeta_  = precondfZeta;
  
  if(verbose_ == true ) 
     cout << "deDSI: precondition factor is set to pf psi: " << precondfPsi_ << " pf theta: " << precondfTheta_<< " pf zeta: " << precondfZeta_ << endl;
}

void DEDSI::setTime( double T , double dtmax, double dtmin, double timestep){

  T_  = T; //  integration time 
  dtmax_ = dtmax;     // integration step size
  dt_    = timestep;
  dtmin_ = dtmin;     // integration step size

  if(verbose_ == true ) cout << " This is set time. T is set to " << T_ << " and dt is set to " << dt_ << endl; 
  return ;
} 

   
void DEDSI::save (const VectorXd & x, const std::string filebase, const std::string outdir, const bool fieldsonly ){
  if(verbose_ == false )
     cout << "!!!!deDSI: save is called" << endl;
  
  
  extractVector(x);      

  np::ndarray u = makeField(x);
  bp::object saveGuess = deModule.attr("save_guess");
  saveGuess(outdir, filebase, u);

  if(rank_ ==0){
     if(Tsearch_)   fout1 << std::setprecision(18) << T_  << endl; 
     if(xrelative_) fout2 << std::setprecision(18) << ax_ << endl; 
     chflow::save(ax_, outdir+"ax_best");
     chflow::save(T_,  outdir+"T_best");
     
  }

}

void DEDSI::saveEigenvec (const VectorXd& x, const string label, const string outdir) {
   cout << "DSI::saveEigenvec() real" << endl;
   extractVector(x); 
   //FIXME: x should be rescaled
   np::ndarray u = makeField(x);
   string filebase = "eigenR"+label;
   bp::object saveGuess = deModule.attr("save_guess");
   saveGuess(outdir, filebase, u);
}

void DEDSI::saveEigenvec (const VectorXd& x1, const VectorXd& x2, const string label1, const string label2, const string outdir) {
   cout << "DSI::saveEigenvec() complex" << endl;
   extractVector(x1); 
   extractVector(x2); 
   //FIXME: x1 & x2 should be rescaled
   np::ndarray u1 = makeField(x1);
   np::ndarray u2 = makeField(x2);
   string filebase1 = "eigenI"+label1;
   string filebase2 = "eigenI"+label2;
   bp::object saveGuess = deModule.attr("save_guess");
   saveGuess(outdir, filebase1, u1);
   saveGuess(outdir, filebase2, u2);
}


Real DEDSI::observable2( VectorXd& x) {
  if(rank_ == 0 && verbose_ == true )
     cout << "deDSI: observable2 (energy) called" << endl;

  np::ndarray u = makeField(x);
   bp::object fEng = deModule.attr("energy");
   bp::object eng = fEng(solver, u);

   Real m = bp::extract<double>(eng);
   return m;
}

Real DEDSI::observable( VectorXd& x) { // FIXME find a better name & position for this function

   if(verbose_ == true ) cout << " deDSI: l2norm observable is called " << endl; 

   np::ndarray u = makeField(x);

   bp::object fl2norm = deModule.attr("l2norm");
   bp::object l2 = fl2norm(solver, u);

   Real m = bp::extract<double>(l2);

   if(rank_ == 0 && verbose_ == true )
      cout << "deDSI: l2norm in observable is " <<scientific << setprecision(17) << m << endl;

   return m;

}

Real DEDSI::DSIL2Norm( const VectorXd& x){
  if(verbose_ == true ) cout << " This is DSIL2Norm " << endl; 
   VectorXd y = x;
   return observable(y);
}

Real DEDSI::DSIL2Norm_X( const VectorXd& x) { // FIXME find a better name & position for this function

   if(verbose_ == true ) cout << " deDSI: l2norm_x is called " << endl; 

   np::ndarray u = vectorXdToNumpyArray(x);

   bp::object fl2norm_x = deModule.attr("l2norm_x");
   bp::object l2_x = fl2norm_x(Nx_, Ny_p_, u);

   Real m = bp::extract<double>(l2_x);

   if(rank_ == 0 && verbose_ == true )
      cout << "deDSI: l2norm_x in observable is " <<scientific << setprecision(17) << m << endl;

   return m;

}


np::ndarray DEDSI::vectorXdToNumpyArray1d( const VectorXd & v){

  int Nx = Nx_;
    
  if(two_fields_ == 1) Nx = 2 * Nx_; 

  bp::tuple shape = bp::make_tuple(Nx);
  np::dtype dtype = np::dtype::get_builtin<double>();
  np::ndarray na  = np::zeros(shape, dtype);

  for ( int nx = 0; nx <Nx; nx++) {
        na[nx] = v(nx);
  }
  
  return na;
}


np::ndarray DEDSI::vectorXdToNumpyArray2d( const VectorXd & v){

  int Nx = Nx_;
  int Ny = Ny_p_;

  if(two_fields_ == 1)   Nx = 2 * Nx_; 
  if(three_fields_ == 1) Nx = 3 * Nx_; 

  bp::tuple shape = bp::make_tuple(Nx , Ny);
  np::dtype dtype = np::dtype::get_builtin<double>();
  np::ndarray na  = np::zeros(shape, dtype);
    

  for( int nx = 0 ; nx < Nx ; nx++) { 
       for( int ny = 0 ; ny < Ny ; ny ++) {
            na[nx][ny] = v(nx * Ny + ny);
       }
  }

  return na;
}


np::ndarray DEDSI::vectorXdToNumpyArray( const VectorXd & v) {

       if(d_ == 1) return vectorXdToNumpyArray1d (v);
  else if(d_ == 2) return vectorXdToNumpyArray2d (v);
  else throw std::invalid_argument( "Only 1d and 2d is implemented" );
}


VectorXd DEDSI::numpyArrayToVectorXd1d(const np::ndarray & na) {

  int Nx = Nx_;

  if(two_fields_ == 1) Nx = 2 * Nx_;


  VectorXd vxd(Nx);

  for( int nx = 0; nx < Nx; nx++){
       vxd(nx) = bp::extract<double>(na[nx]);
  } 

  return vxd;
}


VectorXd DEDSI::numpyArrayToVectorXd2d(const np::ndarray & na) {

  int Nx = Nx_;
  int Ny = Ny_p_;

  if(two_fields_ == 1)   Nx = 2 * Nx_;
  if(three_fields_ == 1) Nx = 3 * Nx_;


  VectorXd vxd(Nx * Ny);

  for ( int nx = 0 ; nx <Nx ; nx++) {
      for( int ny = 0 ; ny < Ny ; ny ++) {
           vxd(nx * Ny + ny) = bp::extract<double>(na[nx][ny]);
      }
  }

  return vxd;
}


VectorXd DEDSI::numpyArrayToVectorXd(const np::ndarray & na) {

       if(d_ == 1) return numpyArrayToVectorXd1d(na); 
  else if(d_ == 2) return numpyArrayToVectorXd2d(na); 
  else throw std::invalid_argument( "Only 1d and 2d is implemented" );
}


VectorXd DEDSI::boostObjectToVectorXd1d(const bp::object & bpo){

  int Nx = Nx_;

  if(two_fields_ == 1)   Nx = 2 * Nx_; 
  if(three_fields_ == 1) Nx = 3 * Nx_; 

  VectorXd vxd(Nx);

  for( int nx = 0 ; nx < Nx ; nx++){
       vxd(nx) = bp::extract<double>(bpo[nx]);
  }

  return vxd;
}


VectorXd DEDSI::boostObjectToVectorXd2d(const bp::object & bpo){

  int Nx = Nx_;
  int Ny = Ny_p_;

  if(  two_fields_ == 1) Nx = 2 * Nx_; 
  if(three_fields_ == 1) Nx = 3 * Nx_; 


  VectorXd vxd(Nx * Ny);

  for( int nx = 0 ; nx < Nx ; nx++) {
       for( int ny = 0 ; ny < Ny ; ny++) {
            vxd(nx * Ny + ny) = bp::extract<double>(bpo[nx][ny]);
       }
  }

   return vxd;
}


VectorXd DEDSI::boostObjectToVectorXd(const bp::object & bpo){

       if(d_ == 1) return boostObjectToVectorXd1d(bpo);
  else if(d_ == 2) return boostObjectToVectorXd2d(bpo);
  else throw std::invalid_argument( "Dimension must be less than 3" );
}


string DEDSI::statsHeader() {

    return "#t\tL2(x)\tEnergy" ;//\tmin(x)\tmax(x)";
}

string DEDSI::stats(const VectorXd& x) {

  stringstream ss;
  string tab ="\t";

  VectorXd x0 = x;
  setprecision(17);
  ss << Ra_ << tab;
  ss << observable(x0) << tab;
  ss << observable2(x0) << tab;
  return ss.str();
}

VectorXd DEDSI::makeVector ( const bp::object & u ){

   if(verbose_ == true ) cout <<" makeVector is called " << endl; 


   VectorXd v = boostObjectToVectorXd(u);
   int uunk = getSize();

   const int Tunk = (  Tsearch_ && rank_ == 0) ? uunk : -1;                       // index for T unknown
   const int xunk = (xrelative_ && rank_ == 0) ? uunk + Tsearch_ : -1;
   int Nunk = (rank_ == 0) ? uunk + Tsearch_ + xrelative_  : uunk;


   if (v.rows() < Nunk)
    v.conservativeResize(Nunk);

   if (rank_ == 0) {
     if (Tsearch_)
       v (Tunk) = T_;
     if (xrelative_)
       v (xunk) = ax_;
   }

  return v;
}


np::ndarray DEDSI::makeField (const VectorXd& x) {

  if(verbose_ == true) cout <<" makeField is called " << endl;

  int uunk = getSize(); // (u); // number of components in x that corresond to u

  VectorXd v(uunk);

  if (rank_ == 0) {
    for(int i = 0; i <uunk; ++i) v(i) = x(i); //TODO in rank 0 I am stripping the Torb and ax from the vector
  }
  else{
   v  = x;
  }


 np::ndarray u = vectorXdToNumpyArray(v);

 return u;

}


void DEDSI::extractVector (const VectorXd& x) {

  if(verbose_ == true ) cout << " extractVector is called " << endl;

  int uunk = getSize(); // (u); // number of components in x that corresond to u

  const int Tunk = uunk + Tsearch_ - 1;
  const int xunk = uunk + Tsearch_ + xrelative_ - 1;
  
  Real ax, T;

  if(verbose_) cout << "!!!!!!!in extract vector rank: " << rank_ << " T: " << T << " ax: " << ax << " uunk " << uunk <<" Tunk " << Tunk << " xunk " << xunk <<endl;

  if (rank_ == 0) {
    T    = Tsearch_   ? x (Tunk) : T_;
    ax   = xrelative_ ? x (xunk) : ax_;
    cout << "in extract vector T is set to " << T << endl;
  }

#ifdef HAVE_MPI
  MPI_Bcast (&ax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&T, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

 if(verbose_ ) cout << "!!!!!!!in extract vector: after bcast rank: " << rank_ << " T: " << T << " ax: " << ax << endl;

 setSigma(sp_, st_, sx_, sy_, ax, T);

}

VectorXd DEDSI::precondition(const VectorXd& x){

 int Nx = Nx_;	
 int Ny = Ny_p_;

 VectorXd v = x;

 if(precondf_ == true) {

   for(int i = 0 ; i < 3 * Nx * Ny ; i++){
            if(i < Nx * Ny )                          v(i) = x(i) * precondfPsi_;
       else if(i > Nx * Ny - 1 && i < 2 * Nx * Ny )   v(i) = x(i) * precondfTheta_;
       else                                           v(i) = x(i) * precondfZeta_;
   } 
 }

 return v;
}

VectorXd DEDSI::unPrecondition(const VectorXd& x){

 int Nx = Nx_;	
 int Ny = Ny_p_;
 
 VectorXd v = x;
 
 if(precondf_ == true) {

   for(int i = 0 ; i < 3 * Nx * Ny ; i++){
            if(i < Nx * Ny )                          v(i) = x(i) / precondfPsi_;
       else if(i > Nx * Ny - 1 && i < 2 * Nx * Ny )   v(i) = x(i) / precondfTheta_;
       else                                           v(i) = x(i) / precondfZeta_;
   } 
 }

 return v;
}

VectorXd DEDSI::eval (const VectorXd& x){


  if(verbose_ == true) cout << " This is eval " << endl; 

  VectorXd gx = G( x );

  int gxsize = gx.size();

  if(rank_ == 0) {
        gx.conservativeResize(x.size());
        for(int i = 0; i < (x.size()- gxsize); i++) gx( gxsize + i) = 0;
  }

  if(verbose_ ) cout << " eval : dsi-l2norm of gx is: " << DSIL2Norm(gx) << " rank: " << rank_ << " size: " << gx.size() << endl;
  silX_ = gx;
  return gx;
} 


VectorXd DEDSI::G (const VectorXd& x) {

   if(verbose_ == true ) cout << " This is G " << endl; 

    extractVector(x); 


    np::ndarray u = makeField(x);
    VectorXd fT_u = f( u , T_);
    

    VectorXd v = numpyArrayToVectorXd(u); 
    VectorXd gv = sigmaOps(fT_u) - v;
    

    return gv; 
}

VectorXd DEDSI::f (const np::ndarray & u , double T) {

  if(verbose_ == false ) cout << " This is f " <<  " T is " << T << " dt_ is " << dt_ << endl; 

  //adjust_for_T(T);

  try {
        integrate = deModule.attr("integrate");
        uT = integrate(solver, T, u, dt_);

      } catch (bp::error_already_set &) {
               handle_pyerror();
      }

  VectorXd fT_u = boostObjectToVectorXd(uT) ;
  

  return fT_u;

}

void DEDSI::setSigma(int sp, int st, int sx, int sy, double ax, double T){

   if(verbose_ == true ) cout << " This is set Sigma " << endl; 

     sp_ = sp;
     st_ = st;
     sx_ = sx;
     sy_ = sy;

     ax_ = ax;
     T_ = T;

   if(verbose_ == true )
      cout << "deDSI: sigma is set to: sp " << sp_<< " st " << st_ << " sx " << sx_ << " sy " << sy_ << " ax  "<< ax_ << " T "<< T_ << endl;

}

bp::object DEDSI::translator(const np::ndarray & a){ //assumes and returns stripped vector 

   if(verbose_ == true ) cout << " This is translator " << endl; 

     np::ndarray fu0 = a;
     
     if(rank_==0 && verbose_ == true )
        cout << "deDSI:translator is called!" << " ax is " << ax_ << endl;

     try {
        translate = deModule.attr("translation");
        sfu = translate(solver, ax_, fu0, 3);
       
      } catch (bp::error_already_set &) {
               handle_pyerror();
      }

     return sfu;
}

bp::object DEDSI::reflector(const np::ndarray & a){ //assumes and returns stripped vector 

   if(verbose_ == true ) cout << " This is reflector " << endl; 

     np::ndarray fu0 = a;

     try {
        reflect = deModule.attr("reflection");
        sfu = reflect(solver, fu0, sp_,st_, sx_, sy_);
       
      } catch (bp::error_already_set &) {
               handle_pyerror();
      }

     return sfu;
}

VectorXd DEDSI::sigmaOps (const VectorXd& fu) {

    bp::object fu1;

    if(verbose_ == true ) cout << " This is sigma ops " << endl; 

    if(verbose_ == true ) cout << "deDSI: sigmaOps is called!" << endl;

    if ( sp_ == 1 && st_ == 1 && sx_ == 1 && sy_ == 1 && xrelative_ == false){

      if(verbose_ == true) cout << " deDSI: identity escape clause is applied "<< endl;
      return fu; // identity escape clause

    }
    else{ 

       np::ndarray fu0 = vectorXdToNumpyArray(fu);
  
       if(ax_ != 0){

         if(verbose_ == false ) cout << "deDSI: dsi will be used to translate the solution " << endl;
     
         fu1 = translator(fu0);

       }
       else if ( (sp_ !=1 || st_ != 1) || (sx_ != 1 || sy_ !=1)){ //FIXME if this should be if or else if

         if(verbose_ == true ) cout << "deDSI: dsi will be used to search for symmetric solutions " << endl;

         fu1 = reflector(fu0);
       }
       
       VectorXd sigma_fu = boostObjectToVectorXd(fu1) ;
       
       cout << " DSIL2Norm(fu) "  << DSIL2Norm(fu) << endl;
       cout << " DSIL2Norm(sigma_fu) "  << DSIL2Norm(sigma_fu) << endl;

       if(verbose_ == true ) cout << "deDSI: size of sigmaOps result is " << sigma_fu.size() << endl;
 
       return sigma_fu;
    } 
}

VectorXd DEDSI::xdiff (const VectorXd& x) {


  if(verbose_ == false ) 
     cout << "deDSI: xdiff is called!" << endl;

  np::ndarray u = makeField(x);

  try {
        xDiff = deModule.attr("xdifferentiate");
        dud_x = xDiff(solver, u ,3);

      } catch (bp::error_already_set &) {
               handle_pyerror();
      }
  
  VectorXd dudx = boostObjectToVectorXd(dud_x) ;


  int dudxsize = dudx.size();

  if(rank_ == 0) { 
        dudx.conservativeResize(x.size());
        for(int i = 0; i < (x.size()-dudxsize);i++) dudx( dudxsize + i) = 0;
  }


   dudx *= 1./chflow::L2Norm (dudx);
   
   return dudx;
//  return gx1;
   
}

VectorXd DEDSI::tdiff (const VectorXd& a, Real epsDt) {
  

  if(verbose_ == false )
     cout << "deDSI: tdiff is called with time " << epsDt << endl;

   np::ndarray u = makeField(a);

   VectorXd v = numpyArrayToVectorXd(u);

   VectorXd dadt = f( u , epsDt ) - v;


   int dadtsize = dadt.size();

   if(rank_ == 0) {dadt.conservativeResize(a.size());
        for(int i = 0; i < (a.size() - dadtsize); i++) dadt(dadtsize + i) = 0;
   }

   dadt *= 1./chflow::L2Norm (dadt);
    
    return dadt;
}

void DEDSI::adjust_for_T (Real T) {

  if(verbose_ == true ) cout << "this is adjust for T " << endl;

  T_ = T;
  if (T < 0) {
    cerr << "TimeStep::adjust_for_T : can't integrate backwards in time.\n"
         << "Exiting." << endl;
    exit (1);
  }
  if (T == 0) {
    dt_ = 0;
    T_ = 0;
  }

  int n = chflow::Greater (iround (T/dt_), 1);
  Real dt = T/n;
  
  while (dt < dtmin_ && n>2 && dt != 0)
    dt = T/--n;
  while (dt > dtmax_ &&  n <= INT_MAX && dt != 0)
    dt = T/++n;
    
    cout << " ADJUST FOR T: " << T << " n " << n << " dt " << dt << endl; 
  

  dt_ = dt;

  if(std::isnan(dt) || std::isnan(T) == true ){ cerr <<" Not a number error!!! dt_: rank:" << rank_ << dt_ << " T_: " << T_  << " n: " << n << endl; exit(1); }
  if(verbose_ == false ) cout << "adjust_for_T: rank "<< rank_ <<" T is set to "<< T_ << " dt is set to : " << dt_ << endl;

}


} // end of namespace
