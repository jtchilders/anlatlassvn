// $Modified: Tue Apr  7 12:51:12 2009 by uwer $
#include "FourMomentum.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <cfloat>
#include <signal.h>

#define USE_SINCOS
#undef MICROSOFT
#undef WITHMESSAGE

#define STEFAN

using namespace std;

int spaversion(){
#ifdef STEFAN
  /* This is the old version which is not defined for momenta parallel to the
     3-direction */
  return(0);
#else
  return(1);
#endif 
}

double FourMomentum::tiny = 1.e-12;

void FourMomentum::printMaple() const {
  cout << name << ":=array(0..3, [" 
       << mom[0] << "," << mom[1] << "," 
       << mom[2] << "," << mom[3] << "]);" << endl;
}

void FourMomentum::printFortran() const {
  cout << name <<"(0) = "<< mom[0] << "d0" << endl;
  cout << name <<"(1) = "<< mom[1] << "d0" << endl;
  cout << name <<"(2) = "<< mom[2] << "d0" << endl;
  cout << name <<"(3) = "<< mom[3] << "d0" << endl;
}  

double FourMomentum::getBeta() const{

  return(sqrt( mom[1]*mom[1] + mom[2]*mom[2] + mom[3]*mom[3] ) / mom[0] );

}

double FourMomentum::kt(const FourMomentum & ref) const {
  /*
   *  Total momentum minus longitudinal gives the transversal:
   *
   */
  double longitudinal = 
    ( ref.mom[1] * mom[1] + ref.mom[2] * mom[2] + ref.mom[3] * mom[3])
    / ref.norm3();
  
  return( sqrt( mom[1]*mom[1] + mom[2]*mom[2] + mom[3]*mom[3]
		 - longitudinal * longitudinal ) );

}

double FourMomentum::kl(const FourMomentum & ref) const {

  return( (ref.mom[1]*mom[1] + ref.mom[2] * mom[2] + ref.mom[3] * mom[3])
    / ref.norm3() );
  
}

double FourMomentum::eta(const FourMomentum & ref) const {
  /*
   * Pseudo rapidity
   */
  double tmp = norm3() + kl(ref);
  double pt = kt(ref);

  return( 0.5 * log( tmp*tmp/(pt*pt)) );
}

double FourMomentum::cos(const FourMomentum & ref) const {
  return( ( ref.mom[1]*mom[1] + ref.mom[2]*mom[2] +ref.mom[3]*mom[3])
	  / norm3() / ref.norm3() );
}

double FourMomentum::phi(const FourMomentum & ref1, const FourMomentum ref2) 
  const {
  FourMomentum hlp = (*this) - dot3p(*this,ref1)*ref1 
    / ( ref1.norm3() * ref1.norm3() );
  return( acos( dot3p(hlp,ref2) / ( hlp.norm3() * ref2.norm3() ) ) );
}

double FourMomentum::y(const FourMomentum & ref) const {
  /*
   * Pseudo rapidity
   */
  double pl = kl(ref);
  double pt = kt(ref);

  return( log( (mom[0]+pl) / sqrt( getmass2() + pt*pt ) ) );

}

double FourMomentum::norm3() const {
  return(sqrt( mom[1]*mom[1] + mom[2]*mom[2] + mom[3]*mom[3] ) ); 
}

void FourMomentum::setFourMomentum(const FourMomentum & p){
  mom[0]=p.mom[0];
  mom[1]=p.mom[1];
  mom[2]=p.mom[2];
  mom[3]=p.mom[3];
}

void FourMomentum::setFourMomentum(const double p[4]){

  mom[0] = p[0]; mom[1] = p[1]; mom[2] = p[2]; mom[3] = p[3];

  if (ext == true) 
    checkmass();

}

void FourMomentum::setFourMomentum(const double p0,const double p1, 
		       const double p2, const double p3){

  mom[0] = p0; mom[1] = p1; mom[2] = p2; mom[3] = p3;

  if (ext == true)
    checkmass();

}

void FourMomentum::boost(const FourMomentum & p, double m){

  m = (m == 0.) ? p.getmass() : m ;
  const double gamma = p.mom[0] / m; 
  const double invp0 = 1.0 / p.mom[0];
 
  const double beta1 = p.mom[1] * invp0;
  const double beta2 = p.mom[2] * invp0;
  const double beta3 = p.mom[3] * invp0;

  const double help = beta1 * mom[1] + beta2 * mom[2] + beta3 * mom[3];
  const double fac = ( gamma / ( 1.0 + gamma ) * help + mom[0] ) * gamma;

  mom[0] = ( mom[0] + help ) * gamma ;
  mom[1] += fac * beta1;
  mom[2] += fac * beta2;
  mom[3] += fac * beta3;

}

bool FourMomentum::checkmass(){

  double tmp1 = mom[0];
  double tmp2 = sqrt( mom[1]*mom[1] + mom[2]*mom[2] + mom[3]*mom[3] );

  if ( mass2 == 0.) {

    double help = fabs(tmp1) + fabs(tmp2);

    if ( fabs(tmp1-tmp2) <= help * tiny  ) 
      return true;

    int prec = cout.precision(15);
    cout << "# Error when setting FourMomenta "<< name <<"."<< endl;
    cout << "# Mass doesn't match stored value." << endl;
    cout << "# Is: m*m = "<< tmp1-tmp2 << " Should: " << mass2 << endl;
    cout.precision(prec);
    return false;
  } 


  if (fabs((tmp1-tmp2)*(tmp1+tmp2) - mass2) < mass2 * tiny  ) 
    return true;

  int prec = cout.precision(15);
  cout << "# Error when setting FourMomenta "<< name <<"."<< endl;
  cout << "# Mass doesn't match stored value." << endl;
  cout << "# Is: m*m = " 
       << (tmp1-tmp2) * (tmp1+tmp2 )  
       << " Should: " << mass2 << endl;
  cout.precision(prec);
  return false;
}

const FourMomentum  operator-(const FourMomentum & p1){
  return( FourMomentum("-(" + p1.name +")", 
		       -p1.mom[0],-p1.mom[1],-p1.mom[2],-p1.mom[3]) );
}

const FourMomentum & operator+(const FourMomentum & p1){
  return(p1);
}

const FourMomentum operator-(const FourMomentum & p1, 
			     const FourMomentum & p2){
  return(FourMomentum( p1.name + "-(" + p2.name +")",
		       p1.mom[0] - p2.mom[0],
		       p1.mom[1] - p2.mom[1],
		       p1.mom[2] - p2.mom[2],
		       p1.mom[3] - p2.mom[3]
		       )
	 );
}

const FourMomentum operator+(const FourMomentum & p1, 
			       const FourMomentum & p2){
  return(FourMomentum( p1.name + "+" + p2.name,
		       p1.mom[0] + p2.mom[0],
		       p1.mom[1] + p2.mom[1],
		       p1.mom[2] + p2.mom[2],
		       p1.mom[3] + p2.mom[3]
		       ));
}

const FourMomentum operator*(const double & s, const FourMomentum & p){
  return(FourMomentum( s*p.mom[0],
		       s*p.mom[1],
		       s*p.mom[2],
		       s*p.mom[3]));
}

const FourMomentum operator*(const FourMomentum & p, const double & s ){
  return(FourMomentum( s*p.mom[0],
		       s*p.mom[1],
		       s*p.mom[2],
		       s*p.mom[3]
		      ));
}

const FourMomentum operator/(const FourMomentum & p,const double & s ){
  return(FourMomentum( p.mom[0]/s,
		       p.mom[1]/s,
		       p.mom[2]/s,
		       p.mom[3]/s
		       ));
}


std::ostream& operator<<(std::ostream& os, const FourMomentum& P){
        return( os << P.name <<"(0) = "<< P.mom[0] << ";\n"
	  << P.name <<"(1) = "<< P.mom[1] << ";\n"
	  << P.name <<"(2) = "<< P.mom[2] << ";\n"
	  << P.name <<"(3) = "<< P.mom[3]<< ";\n");
}


void cross3p(const FourMomentum& p1, const FourMomentum& p2,
	     FourMomentum& p3){

  p3.mom[1] = p1.mom[2]*p2.mom[3] - p1.mom[3]*p2.mom[2];
  p3.mom[2] = p1.mom[3]*p2.mom[1] - p1.mom[1]*p2.mom[3];
  p3.mom[3] = p1.mom[1]*p2.mom[2] - p1.mom[2]*p2.mom[1];

}

double dot3p(const FourMomentum & p1, const FourMomentum & p2){
  /*
   * Takes the space components and evaluates the dotproduct
   */
  return( p1.mom[1]*p2.mom[1] + p1.mom[2]*p2.mom[2] + p1.mom[3]*p2.mom[3] );
}

double dotp(const FourMomentum & p1, const FourMomentum & p2) {

  if ( &p1 == &p2 )
    return( p1.getmass2() );

  double t2 =  p1.mom[0] * p2.mom[0];
  double x2 =  p1.mom[1] * p2.mom[1]  
    + p1.mom[2]*p2.mom[2] 
    + p1.mom[3]*p2.mom[3];
  
  double tmp = t2 - x2;

  if (fabs(tmp/t2) < DBL_EPSILON)
    return 0.0;

  return(tmp);

}

void rotate3(double angle, FourMomentum & n, FourMomentum & x) {
  /*
   * Brute force routine to rotate a vector arround n.
   * Not optimized for anything --- it should just work.
   */
  const double sinphi = sin(angle);
  const double cosphi = cos(angle);
  static double e[3][3][3];
  static bool initialized = false;

  if ( !initialized ) {
    for ( int i = 0; i < 3; i++)
      for ( int j = 0; j < 3; j++)
	for ( int k = 0; k < 3; k++)
	  e[i][j][k] = ( i == j || j == k || k == i ) ? 0 :
	    ( (i-j)*(j-k)*(k-i) > 0 ) ? 1 : -1;
    initialized = true;
  }

  double R[4][4];

  for (unsigned int i = 1; i < 4 ; i++){
    for (unsigned int j = 1; j < 4 ; j++){
      R[i][j] = n.mom[i]*n.mom[j] * ( 1.0 - cosphi )
	- ( + n.mom[1] * e[0][i-1][j-1] 
	    + n.mom[2] * e[1][i-1][j-1] 
	    + n.mom[3] * e[2][i-1][j-1] ) * sinphi; 
    }
  }

  for (unsigned int i = 1; i < 4 ; i++){
    R[i][i] += cosphi;
  }

  for (unsigned int i = 1; i < 4 ; i++){
    for (unsigned int j = 1; j < 4 ; j++){
      
    }
  }
  
  double tmp[4];

  for (unsigned int i = 1; i < 4; i++) {
    tmp[i] = 0;
    for (unsigned int j = 1; j < 4; j++) {
      tmp[i] += R[i][j] * x.mom[j];
    }
  }

  for (unsigned int i = 1; i < 4; i++) {
    x.mom[i] = tmp[i];
  }

} 

FourMomentum boost(const FourMomentum & p, const FourMomentum & x){
  FourMomentum res(x);
  res.boost(p);
  return(res);
}

std::complex<double> spa(const FourMomentum & p, 
			 const FourMomentum & q){
  /*
   *    Formualae are taken from  Stefans thesis
   */

#ifdef STEFAN
  double pplus = p.mom[0] + p.mom[3];
  double qplus = q.mom[0] + q.mom[3];
#else
  double pplus = p.mom[0] + p.mom[1];
  double qplus = q.mom[0] + q.mom[1];
#endif  
  /*
   * We use:
   *
   * a = p.mom[1]*qplus + complex<double>(0.0,1.0) * qplus * p.mom[2];
   * b = q.mom[1]*pplus + complex<double>(0.0,1.0) * pplus * q.mom[2];
   *
   * spa = (0.d0,-1.d0)**n / dsqrt(dabs(pplus*qplus)) * (a - b)
   */

  /*
   * complex<double> a_minus_b( p.mom[1]*qplus - q.mom[1]*pplus,
   *			     qplus * p.mom[2] - pplus * q.mom[2]);
   *  complex<double> help = 1.0 / sqrt( fabs( pplus*qplus ) ) * ( a_minus_b );
   */

  //cout<<pplus<<" "<<qplus<<endl;
  double prefac = 1.0 / sqrt( fabs( pplus*qplus ) );

#ifdef STEFAN
  double x = ( qplus * p.mom[1] - pplus * q.mom[1] ) * prefac;
  double y = ( qplus * p.mom[2] - pplus * q.mom[2] ) * prefac;
#else
  double x = ( qplus * p.mom[2] - pplus * q.mom[2] ) * prefac;
  double y = ( qplus * p.mom[3] - pplus * q.mom[3] ) * prefac;
#endif
  complex<double> result[] = { complex<double>(x,y),
			       complex<double>(y,-x),
			       complex<double>(-x,-y)};

  int n = 0;
  if ( p.mom[0] < 0.0 ) n=1;
  if ( q.mom[0] < 0.0 ) n++;

  /*
   * NOTE: to avoid branchings we may replace the above code
   *       by:
   *
   * #define DOUBLESIGN2INT(_X_) (( * ( (int*)&_X_ + 1)  & 0x80000000 )>> 31 )
   * int n = DOUBLESIGN2INT(p.mom[0]) + DOUBLESIGN2INT(q.mom[0]);
   *
   */
  return( result[n] );

}


std::complex<double> spb(const FourMomentum & p1, 
			 const FourMomentum & p2){
  return( p1.mom[0]*p2.mom[0] < 0.0 ? conj(spa(p1,p2)) 
	  : conj(spa(p2,p1) ) );
}

void SplitMomenta(const FourMomentum & kq, FourMomentum & q1, 
		  FourMomentum & q2){
  /* 
   * Routine to splitt massive momentum kq into two massless momenta q1,q2
   *
   * We use the helicity basis to define q1, q2. 
   */
  double m = kq.getmass();
  double n = m/kq.norm3();

  FourMomentum sq(0., kq.mom[1]*n, kq.mom[2]*n, kq.mom[3]*n);

  sq.boost(kq);

  for( int i = 0; i < 4; i++ ){
    q1.mom[i] = 0.5 * ( kq.mom[i] + sq.mom[i] ) ;
    q2.mom[i] = 0.5 * ( kq.mom[i] - sq.mom[i] ) ;
  }


}
 

void evalqr(const FourMomentum & kq, const FourMomentum & kqb, 
            FourMomentum & q1, FourMomentum & q2, 
            FourMomentum & r1, FourMomentum & r2){

  int err = 0;

  double m = kq.getmass();

  double nq  = m/kq.norm3();
  double nqb = m/kqb.norm3();

  /*
   *  ------ Spin in the rest frames ------------------------
   * We normalize the spin vector to m^2 to save a few multiplications.
   */
  FourMomentum sq(0.0,  kq.mom[1]*nq, kq.mom[2]*nq, kq.mom[3]*nq );
  FourMomentum sqb(0.0, kqb.mom[1]*nqb, kqb.mom[2]*nqb, kqb.mom[3]*nqb );

  /*
   *  ------ Calculate the spins in the cms: ---------------             
   */
  sq.boost(kq);
  sqb.boost(kqb); 
  
  /*            
   * q1 = ( kq + m * sq ) / 2.0;
   * q2 = ( kq - m * sq ) / 2.0;
   *
   *  ------ Note: the relation is reversed for the anti-quark:
   *
   * r1 = ( kqb - m * sqb ) / 2.0;
   * r2 = ( kqb + m * sqb ) / 2.0;
   * 
   * The factor m is already absorbed into the definition of s.
   */
  for( int i = 0; i < 4; i++ ){
    q1.mom[i] = 0.5 * ( kq.mom[i] + sq.mom[i] );
    q2.mom[i] = 0.5 * ( kq.mom[i] - sq.mom[i] );
    r1.mom[i] = 0.5 * ( kqb.mom[i] - sqb.mom[i] );
    r2.mom[i] = 0.5 * ( kqb.mom[i] + sqb.mom[i] );
  }

}


static inline double lambda(const double &a, const double &b, const double &c){
  return a*a + b*b + c*c - 2.0*a*b - 2.0*b*c - 2.0*a*c;
} 

int event(double x[], vector<FourMomentum> & kout, double & jacobi) {
  
  FourMomentum ktmp[10];
  
  int n = kout.size()-2;

  for (int i = 0; i < n; i++){
    ktmp[i] = kout[i+2];
  }

  double s = dotp(kout[0],kout[0]) + 2.0 * dotp(kout[0],kout[1])
    + dotp(kout[1],kout[1]);

  int retval = eventn(n,s,x,ktmp,jacobi);

  for (int i = 0; i < n; i++){
    kout[i+2] = ktmp[i];
  }
  
  return(retval);
  
}

int eventn(const int n, const double s, const double x[], 
           FourMomentum kout[], double & jacobi){
  
  /*   Routine to generate a n-parton event with center-of-mass
   *   energy squared given by s and the parton masses set through
   *   the 4-momenta in kout[].
   *   The momentum configuration is returned in kout[], the jacobian for
   *   the transformation is returned through jacobi.
   *   x should contain 3*n-4 random numbers
   */
  
  static const double TwoPi = 2.0 * M_PI;
  int err = 0;
  const double sqrts = sqrt(s);

#define MICROSOFT
#ifdef MICROSOFT
  /*
   * The Microsoft compiler doesn't allow to use n in the array definition.
   */
#define MAXN 8

  if ( n > MAXN ) {
    std::cout << "eventn: n too large, redefine MAXN and compile a "
	      << "new version..." << std::endl;
    exit(1);
  }

#else
#define MAXN n
#endif

  //  double pset[MAXN][4] ;
  //  double m[MAXN],MM[MAXN], MMQ[MAXN], mu[MAXN], P[MAXN];

  struct {
    double pset[4];
    double m;
    double MM;
    double MMQ;
    double mu;
    double P;
  } pinfo[MAXN];

  int irn = 0;

  //#define XRANDOM x[irn++]
#define XRANDOM  *(x++)

#ifdef WITHMESSAGE

  static bool initialized = false;
  
  if ( !initialized ) {
    std::cout << "#  -------------------------------------" << std::endl
              << "#  Simple event generator" << std::endl
              << "#  written by P.Uwer " << std::endl
              << "#" << std::endl 
              << "#  [0,1]^(3n-4) --> (p1,...,pn)         " << std::endl
              << "#  based on Byckling&Kajantie Chap. 10   " << std::endl
              << "#  -------------------------------------" << std::endl;
    initialized = true;
  };

#endif

  /*
  for ( unsigned int i = 0; i < n; i++) {
    mu[i] = 0.0;
     m[i] = kout[i].getmass();
    for ( unsigned int j = 0; j <= i; j++) {
      mu[i] += m[j];
    };
  };
  */
  double msum=0;
  for ( unsigned int i = 0; i < n; i++) {
    pinfo[i].m = kout[i].getmass();
    msum += pinfo[i].m;
    pinfo[i].mu = msum;
  };
  
  /*
   *     Set the invariant masses MMQ(i) = (p_0+...+p_i)^2
   */     
  pinfo[n-1].MM  = sqrts; 
  pinfo[n-1].MMQ = s;
  pinfo[0].MM    = pinfo[0].m;
  pinfo[0].MMQ   = pinfo[0].MM * pinfo[0].MM;

  for (unsigned int i = n-2; i > 0; i--) {
    pinfo[i].MM = pinfo[i].mu + XRANDOM * ( pinfo[i+1].MM - pinfo[i+1].mu );
    pinfo[i].MMQ = pinfo[i].MM * pinfo[i].MM;
  }

  /*
   *     Calculate the three momenta and the jacobian
   */
  pinfo[1].P = pinfo[1].m == 0.0 ? ( pinfo[1].MM - pinfo[0].MM ) * 
    ( pinfo[1].MM + pinfo[0].MM ) / ( 2.0 * pinfo[1].MM ) :
      sqrt( lambda(pinfo[1].MMQ, pinfo[0].MMQ, pinfo[1].m*pinfo[1].m ) ) 
    / ( 2.0 * pinfo[1].MM );

  jacobi = pow(TwoPi, n-1) / ( 2.0 * sqrts ) * pinfo[1].P;


  for( unsigned int i = 2; i < n; i++ ) {
    pinfo[i].P = 
      ( pinfo[i].m == 0.0 ? 
	( pinfo[i].MM - pinfo[i-1].MM ) * ( pinfo[i].MM + pinfo[i-1].MM ) 
	/ ( 2.0 * pinfo[i].MM ) :
	sqrt( lambda(pinfo[i].MMQ, pinfo[i-1].MMQ, pinfo[i].m*pinfo[i].m ) ) 
	/ ( 2.0 * pinfo[i].MM ) );
    jacobi *= pinfo[i].P * ( pinfo[i].MM - pinfo[i].mu );
 }

  /*
   *     Finally reconstruct the event in the overall CMS:
   *
   *     Calculate first the 4-momenta pset(rho,i) in the restframe of
   *     p_1+...+p_i:
   */
  for(unsigned int i = 1; i < n; i++) {
#ifdef USE_SINCOS
    double sinphi_i, cosphi_i;
    sincos(TwoPi * XRANDOM, &sinphi_i, &cosphi_i);
#else
    /*
     *  Note: For phi we have to calculate the cosine explicitly!!!!
     *  dsqrt(1-sin^2) would only give cos() > 0
     */
    double phi_i = TwoPi * XRANDOM;
    double sinphi_i = sin( phi_i );
    double cosphi_i = cos( phi_i );
#endif
    double tmp = XRANDOM;
    double costheta_i = 2 * tmp - 1.0;
    double sintheta_i = 2 * sqrt( tmp * ( 1. - tmp ) );
    //double sintheta_i = sqrt(1.0 - costheta_i*costheta_i);

    pinfo[i].pset[1] = pinfo[i].P * cosphi_i*sintheta_i;
    pinfo[i].pset[2] = pinfo[i].P * sinphi_i*sintheta_i;
    pinfo[i].pset[3] = pinfo[i].P * costheta_i;
    pinfo[i].pset[0] = pinfo[i].m == 0. ? pinfo[i].P : 
      sqrt( pinfo[i].m*pinfo[i].m + pinfo[i].P*pinfo[i].P );

    kout[i].setFourMomentum(pinfo[i].pset);
  }

  pinfo[0].pset[0] = pinfo[1].MM - pinfo[1].pset[0];
  pinfo[0].pset[1] = - pinfo[1].pset[1];
  pinfo[0].pset[2] = - pinfo[1].pset[2];
  pinfo[0].pset[3] = - pinfo[1].pset[3];

  kout[0].setFourMomentum(pinfo[0].pset);


  FourMomentum bhelp;
  for (unsigned int i = 2; i < n; i++) {
    bhelp.setFourMomentum( pinfo[i].MM-pinfo[i].pset[0], 
			   -pinfo[i].pset[1], 
			   -pinfo[i].pset[2], 
			   -pinfo[i].pset[3]);
    for (unsigned int j = 0; j < i; j++){
      kout[j].boost(bhelp, pinfo[i-1].MM);
    }
  }  
  return err;
}

int event(double s, double x[2], FourMomentum & p1, FourMomentum & p2,
      double & jacobi ){
  FourMomentum p[]={p1,p2};
  int err;



  err= eventn(2, s, x, p, jacobi);
  p1 = p[0];
  p2 = p[1];
  return(err);
}

int event(double s, double x[5], FourMomentum & p1, FourMomentum & p2,
      FourMomentum & p3,
      double & jacobi ){
  FourMomentum p[]={p1,p2,p3};
  int err;



  err = eventn(3, s, x, p, jacobi);
  p1 = p[0];
  p2 = p[1];
  p3 = p[2];
  return(err);
}

int event(double s, double x[8], FourMomentum & p1, FourMomentum & p2,
      FourMomentum & p3,FourMomentum & p4,
      double & jacobi ){
  FourMomentum p[]={p1,p2,p3,p4};
  int err;



  err= eventn(4, s, x, p, jacobi);
  p1 = p[0];
  p2 = p[1];
  p3 = p[2];
  p4 = p[3];
  return(err);
}

int event(double s, double x[11], FourMomentum & p1, FourMomentum & p2,
      FourMomentum & p3,FourMomentum & p4,FourMomentum & p5,
      double & jacobi ){
  FourMomentum p[]={p1,p2,p3,p4,p5};
  int err;



  err= eventn(5, s, x, p, jacobi);
  p1 = p[0];
  p2 = p[1];
  p3 = p[2];
  p4 = p[3];
  p5 = p[4];
  return(err);
}

#ifdef RAMBOEVENT
extern "C" void phpoint_(const int & n, const double& et, const double xm[], 
	      const double rn[100][4], 
	      double p[100][4], double & jacoobi);

void RamboEvent(double s, double x[8], FourMomentum & p1, FourMomentum & p2,
	   FourMomentum & p3,FourMomentum & p4,
	   double & jacobi ){
  double et = sqrt(s);
  double xm[4];
  double p[100][4];
  
  xm[0] = p1.getmass();
  xm[1] = p2.getmass();
  xm[2] = p3.getmass();
  xm[3] = p4.getmass();

  phpoint_(4,et,xm, (double (*)[4])x, p, jacobi);

  p1.setFourMomentum(p[0][3],p[0][0],p[0][1],p[0][2]);
  p2.setFourMomentum(p[1][3],p[1][0],p[1][1],p[1][2]);
  p3.setFourMomentum(p[2][3],p[2][0],p[2][1],p[2][2]);
  p4.setFourMomentum(p[3][3],p[3][0],p[3][1],p[3][2]);

}

#endif
