// $Modified: Thu Nov 23 00:48:27 2006 by puwer $
#include "StandardModelParameters.h"
#include "Reduction.h"
#include "GGReduction.h"
#include "ScalarInt.h"
#include "FF.h"
#include <complex>
#include <cmath>
/*
 * TODO: check the setting of xd0
 */

#define PRINT(x_) cout << #x_ << " " << x_ << endl;

using namespace std;


const bool OnlyUVSing = false;


const double pi2 = M_PI * M_PI;


static StandardModelParameters& parms = StandardModelParameters::instance();

extern "C" {
complex<double> b0check_(const double & divergence, const double & mudim2, 
	       const double & p, const double & m1, const double & m2);
complex<double> db0check_(const double & p, const double & m1, 
			  const double & m2);
}
/***
 *** Interface to the Li2 fortran routine. Note that for x>1 only the
 *** real part is calculated.
 ***/
extern "C" {
  double li2_(const double & x);
}

/*
 * We consider a momentum squared as massless, if p^2 < tiny.
 * In principle we could ask directly for 0.0, because dotp is written
 * in such a way that if a cancellation of x digits is observed the 
 * the momentum squared is set to zero...
 */
const double tiny=1.e-7;
const double deltam2=1.e-7;


static bool initialized = false;

//// addedy by SA to check the errors
void finishIntegralLibrary(){
  cout << "#\n";
  cout << "# Cleaning the library for the scalar integrals...\n";
  cout << "#\n";
  cout.flush(); 
  ffexi_();
}



void initIntegralLibrary(){
  cout << "#\n";
  cout << "# Init the library for the scalar integrals...\n";
  cout << "#\n";
  cout.flush(); 
  ffini_();
  // Set verbose warning output
  extern bool FFTestFlag;
  if(FFTestFlag) {
    cout << "#\n";
    cout << "# Turning on ltest in ff library...\n";
    cout << "#\n";
    ffchangeflags_();
  }
  //
  if (OnlyUVSing)
    cout << "# Calculate only contrib. propotional to 1/e\n";
  /*
   * We use FF only for IR finite integrals, for all other integrals
   * it should crash..., so we set the regulator mass to zero.
   */
  ffcut_.delta = 0.0;
  initialized = true;
}

Dinput fixDm0(const Dinput & in){

  /* 
   * If there is a zero mass we rotate it to positon 0.
   *
   * The masses are set explicitly so we can ask for m==0...
   */
  if (in.m0q == 0.0) {
    return(in);
  }

  if (in.m1q == 0.0) {
    return(Dinput(in.p2q,in.p3q,in.p6q,in.p1q,in.p7q,in.p5q,
		  in.m1q,in.m2q,in.m3q,in.m0q));
    //    return(Dinput(in.p2,in.p3,-in.p6,in.m1q,in.m2q,in.m3q,in.m0q));
  }

  if (in.m2q == 0.0) {
    return(Dinput(in.p3q,in.p6q,in.p1q,in.p2q,in.p5q,in.p7q,
		  in.m2q,in.m3q,in.m0q,in.m1q));
    //    return(Dinput(in.p3,-in.p6,in.p1,in.m2q,in.m3q,in.m0q,in.m1q));
  }

  if (in.m3q == 0.0) {
    return(Dinput(in.p6q,in.p1q,in.p2q,in.p3q,in.p7q,in.p5q,
		  in.m3q,in.m0q,in.m1q,in.m2q));
    //    return(Dinput(-in.p6,in.p1,in.p2,in.m3q,in.m0q,in.m1q,in.m2q));
  }

  cout << "fixDm0: Something is wrong: tried to move zero mass to pos. 1, \n "
       << " but there is no zero mass " << endl;

  exit(1);

}

Dinput fixDm1(const Dinput & in){
  /* 
   * If there is a second mass which is zero we rotate it to positon 1.
   *
   * Note that we assume here m0q == 0 already.
   *
   * The masses are set explicitly so we can ask for m==0...
   */  
  if (in.m1q == 0.0) {
    return(in);
  }

  if (in.m2q == 0.0) {
    return(Dinput(in.p5q,in.p2q,in.p7q,in.p6q,in.p1q,in.p3q,
		  in.m0q,in.m2q,in.m1q,in.m3q));
    //    return(Dinput(in.p5,-in.p2,in.p7,in.m0q,in.m2q,in.m1q,in.m3q));
  }

  if (in.m3q == 0.0) {
    return(Dinput(in.p6q,in.p3q,in.p2q,in.p1q,in.p5q,in.p7q,
		  in.m0q,in.m3q,in.m2q,in.m1q));
    //    return(Dinput(in.p6,-in.p3,-in.p2,in.m0q,in.m3q,in.m2q,in.m1q));
  }

  cout << "fixDm1: Something is wrong: tried to move zero mass to pos. 1, \n "
       << " but there is no zero mass " << endl;

  exit(1);
}

Dinput fixDm123(const Dinput & in){
  /*
   * If there are 3 masses which are zero, we rotate the non-zero mass 
   * to position 3.
   *
   * Again we can check for mq==0.0 because they are set explicitly.
   *
   */
  if (in.m3q != 0.0) {
    return(in);
  }

  if (in.m0q != 0.0) {
    return(Dinput(in.p3q,in.p2q,in.p1q,in.p6q,in.p7q,in.p5q,
		  in.m3q,in.m2q,in.m1q,in.m0q));
    //    return(Dinput(-in.p3,-in.p2,-in.p1,in.m3q,in.m2q,in.m1q,in.m0q));
  }

  if (in.m1q != 0.0) {
    return(Dinput(in.p6q,in.p3q,in.p2q,in.p1q,in.p5q,in.p7q,
	     in.m0q,in.m3q,in.m2q,in.m1q));
    //    return(Dinput(in.p6,-in.p3,-in.p2,in.m0q,in.m3q,in.m2q,in.m1q));
  }

  if (in.m2q != 0.0) {
    return(Dinput(in.p1q,in.p7q,in.p3q,in.p5q,in.p6q,in.p2q,
	     in.m0q,in.m1q,in.m3q,in.m2q));
    //    return(Dinput(in.p1,in.p7,-in.p3,in.m0q,in.m1q,in.m3q,in.m2q));
  }

  cout << "fixDm123: Something is wrong: tried to move massive mass "
       << "to pos. 3, \n "
       << " but there is no non-zero mass " << endl;
  exit(1);
}


Dinput fixM0(const Dinput & in){
  /*
   * Rotate the massless legs to position 1,2,3.
   */

  if  ( ( fabs(in.p1q) < tiny ) && ( fabs(in.p2q) < tiny ) 
	&& ( fabs(in.p3q) < tiny ) )
    return(in);

  if  ( ( fabs(in.p1q) < tiny ) && ( fabs(in.p3q) < tiny ) 
	&& ( fabs(in.p5q) < tiny ) )
    return(Dinput(in.p1q,in.p5q,in.p3q,in.p7q,in.p2q,in.p6q,
		  in.m1q,in.m0q,in.m2q,in.m3q));
    //    return(Dinput(-in.p1,in.p5,in.p3,in.m1q,in.m0q,in.m2q,in.m3q));

  // added: 26/09/2006
  if  ( ( fabs(in.p1q) < tiny ) && ( fabs(in.p3q) < tiny ) 
	&& ( fabs(in.p6q) < tiny ) )
    //    return(Dinput(in.p3q,in.p6q,in.p1q,in.p2q,in.p5q,in.p7q,
    //	  in.m2q,in.m3q,in.m0q,in.m1q));
    return(Dinput(in.p1q,in.p6q,in.p3q,in.p2q,in.p7q,in.p5q,
		  in.m1q,in.m0q,in.m3q,in.m2q));

  if  ( ( fabs(in.p1q) < tiny ) && ( fabs(in.p5q) < tiny ) 
	&& ( fabs(in.p7q) < tiny ) )
    return(Dinput(in.p5q,in.p1q,in.p7q,in.p3q,in.p2q,in.p6q,
		  in.m2q,in.m0q,in.m1q,in.m3q));
    //    return(Dinput(-in.p5,in.p1,in.p7,in.m2q,in.m0q,in.m1q,in.m3q));

  cout << "fixM0: something went wrong while adjusting M0, exit...\n";
  in.printmaple();
  exit(1);
}

Dinput fixM1a1(const Dinput & in){
  /*
   * Rotate the integral such that p1^2 = p2^2 = 0,
   * p3^2=m^2 is achieved by fixM1a2
   */
  if (( fabs(in.p1q) < tiny)  && ( fabs(in.p2q) < tiny) )
    return(in);


  if ( ( fabs(in.p1q) < tiny ) && ( fabs(in.p5q) < tiny ))
    return(Dinput(in.p1q,in.p5q,in.p3q,in.p7q,in.p2q,in.p6q,
		  in.m1q,in.m0q,in.m2q,in.m3q)); 
    //    return(Dinput(-in.p1,in.p5,in.p3,in.m1q,in.m0q,in.m2q,in.m3q)); 

  if ( ( fabs(in.p2q) < tiny ) && ( fabs(in.p5q) < tiny ) )
    return(Dinput(in.p5q,in.p2q,in.p7q,in.p6q,in.p1q,in.p3q,
	     in.m0q,in.m2q,in.m1q,in.m3q));
    //    return(Dinput(in.p5,-in.p2,in.p7,in.m0q,in.m2q,in.m1q,in.m3q));
  
  cout << "fixM1a1: something went wrong while adjusting M1a, exit...\n";
  in.printmaple();
  exit(1);

}

Dinput fixM1a2(const Dinput & in){
  /*
   * p1q=0, p2q=0, m0=m1=m2=0, make sure that p3q=m3q
   */
  if ( fabs( in.p3q - in.m3q ) < deltam2 )
    return(in);

  if ( fabs( in.p6q - in.m3q ) < deltam2 )
    return(Dinput(in.p2q,in.p1q,in.p6q,in.p3q,in.p5q,in.p7q,
		  in.m2q,in.m1q,in.m0q,in.m3q));
    //    return(Dinput(-in.p2,-in.p1,in.p6,in.m2q,in.m1q,in.m0q,in.m3q));

  cout << "fixM1a2: something went wrong while adjusting M3a, exit...\n";
  in.printmaple();
  exit(1);
}

Dinput fixM1b1(const Dinput & in){
  /*
   * Rotate the massless leg to position 1.
   */
  if ( fabs(in.p1q) < tiny )
    return(in);

  if ( fabs(in.p2q) < tiny )
    return(Dinput(in.p2q,in.p5q,in.p6q,in.p7q,in.p1q,in.p3q,
		  in.m1q,in.m2q,in.m0q,in.m3q)); 
    //    return(Dinput(in.p2,-in.p5,in.p6,in.m1q,in.m2q,in.m0q,in.m3q)); 

  if ( fabs(in.p5q) < tiny )
    return(Dinput(in.p5q,in.p2q,in.p7q,in.p6q,in.p1q,in.p3q,
		  in.m0q,in.m2q,in.m1q,in.m3q)); 
    //    return(Dinput(in.p5,-in.p2,in.p7,in.m0q,in.m2q,in.m1q,in.m3q)); 

  cout << "fixM1b1: something went wrong while adjusting M1b, exit...\n";
  in.printmaple();
  exit(1);
}

Dinput fixM1b2(const Dinput & in){
  /*
   * Second stage to fix M1b, m0=m1=m2=0, p1q = 0, fix now p3^2=p4^2=m^2
   */
  if ( (fabs(in.p3q-in.m3q) < deltam2 ) && (fabs(in.p6q-in.m3q) < deltam2 ) )
    return(in);

  if ( (fabs(in.p3q-in.m3q) < deltam2 ) && (fabs(in.p7q-in.m3q) < deltam2 ) )
    return(Dinput(in.p1q,in.p5q,in.p3q,in.p7q,in.p2q,in.p6q,
		  in.m1q,in.m0q,in.m2q,in.m3q));
    //return(Dinput(-in.p1,in.p5,in.p3,in.m1q,in.m0q,in.m2q,in.m3q));


  cout << "fixM1b2: something went wrong while adjusting M3a, exit...\n";
  in.printmaple();
  exit(1);

}


Dinput fixM3a1(const Dinput & in){
  if (in.m0q != 0.0) {
    cout << "fixM3a: m0 should be already zero, exit program";
    exit(1);
  } 

  /*
   * We rotate first the massless external leg to pos 3.
   */
  if ( fabs(in.p3q) < tiny ) 
    return(in);

  if ( fabs(in.p2q) < tiny ) {
    return(Dinput(in.p6q,in.p7q,in.p2q,in.p5q,in.p1q,in.p3q,
		  in.m0q,in.m3q,in.m1q,in.m2q));
    //return(Dinput(in.p6,-in.p7,in.p2,in.m0q,in.m3q,in.m1q,in.m2q));
  }

  if ( fabs(in.p7q) < tiny ) {
    return(Dinput(in.p5q,in.p2q,in.p7q,in.p6q,in.p1q,in.p3q,
		  in.m0q,in.m2q,in.m1q,in.m3q));
    //return(Dinput(in.p5,-in.p2,in.p7,in.m0q,in.m2q,in.m1q,in.m3q));
  }
 
  cout << "fixM3a1: Something is wrong: tried to move zero"
       << " mass leg to pos. 3, \n "
       << " but there is no zero mass leg" << endl;
  in.printmaple();
  exit(1);
}

Dinput fixM3a2(const Dinput & in){
  /*
   * m0 = 0 and p3q = 0 already, now we rotate the two massive (on-shell)
   * legs to pos 1 and 4.
   */
  if ( ( fabs(in.m1q-in.p1q) < deltam2 ) && ( fabs(in.m1q-in.p6q) < deltam2 ) ){
    // Nothing to do:
    return(in);
  }

  if ( (fabs(in.m1q-in.p1q) < deltam2 ) && ( fabs(in.m1q-in.p5q) < deltam2  ) ){
    return(Dinput(in.p1q,in.p7q,in.p3q,in.p5q,in.p6q,in.p2q,
		  in.m0q,in.m1q,in.m3q,in.m2q));
    //return(Dinput(in.p1,in.p7,-in.p3,in.m0q,in.m1q,in.m3q,in.m2q));
  }

  cout << "fixM3a2: Something is wrong: tried to move massive"
       << " legs to pos. 1 and 4 \n ";
  cout << "# of massless legs " << 6-in.getExtMassCount() << endl;
  cout << in.p3q << endl;
  in.printmaple();
  exit(1);

}

Dinput fixM3b(const Dinput & in){
  if (in.m0q != 0.0) {
    cout << "fixM3b: m0 should be already zero, exit program";
    exit(1);
  } 

  if (( fabs(in.p2q) < tiny ) && ( fabs(in.p3q) < tiny ) )
    if (fabs(in.p1q-in.m2q) < deltam2)
      return(in);
    else
      return(Dinput(in.p6q,in.p3q,in.p2q,in.p1q,in.p5q,in.p7q,
		    in.m0q,in.m3q,in.m2q,in.m1q));
      //return(Dinput(in.p6,-in.p3,-in.p2,in.m0q,in.m3q,in.m2q,in.m1q));

  if (( fabs(in.p2q) < tiny ) && ( fabs(in.p7q) < tiny )) {
    if (fabs(in.p5q-in.m2q) < deltam2 )
      return(Dinput(in.p5q,in.p2q,in.p7q,in.p6q,in.p1q,in.p3q,
		    in.m0q,in.m2q,in.m1q,in.m3q));
    //return(Dinput(in.p5,-in.p2,in.p7,in.m0q,in.m2q,in.m1q,in.m3q));
    else
      return(Dinput(in.p6q,in.p7q,in.p2q,in.p5q,in.p1q,in.p3q,
		    in.m0q,in.m3q,in.m1q,in.m2q));
    //return(Dinput(in.p6,-in.p7,in.p2,in.m0q,in.m3q,in.m1q,in.m2q));
  }

  if ( ( fabs(in.p3q) < tiny ) && ( fabs(in.p7q) < tiny)) {
    if ( fabs( in.p5q-in.m2q ) < deltam2 )
      return( Dinput(in.p5q,in.p3q,in.p7q,in.p1q,in.p6q,in.p2q,
		     in.m0q,in.m2q,in.m3q,in.m1q));
      //return(Dinput(in.p5,in.p3,-in.p7,in.m0q,in.m2q,in.m3q,in.m1q));
    else
      return(Dinput(in.p1q,in.p7q,in.p3q,in.p5q,in.p6q,in.p2q,
		    in.m0q,in.m1q,in.m3q,in.m2q));
      //return(Dinput(in.p1,in.p7,-in.p3,in.m0q,in.m1q,in.m3q,in.m2q));
  }

  cout << "fixM3b: Something went wrong when fixing M3b\n";
  in.printmaple();
  exit(1);

}


Dinput fixM2a(const Dinput & in){

  if ( (in.m0q != 0.0) || (in.m1q != 0.0) ){
    cout << "fixM2a: m0,m1 should be already zero, exit program";
    exit(1);
  } 

  if (( fabs(in.p1q) < tiny ) && ( fabs(in.p3q) < tiny ) )
    if ( fabs( in.p2q - in.m3q ) < deltam2 )
      return(in);
    else {
      if (fabs( in.p7q - in.m3q ) < deltam2 )
	return(Dinput(in.p1q,in.p7q,in.p3q,in.p5q,in.p6q,in.p2q,
		 in.m0q,in.m1q,in.m3q,in.m2q));
	//return(Dinput(in.p1,in.p7,-in.p3,in.m0q,in.m1q,in.m3q,in.m2q));
      if (fabs( in.p5q - in.m3q ) < deltam2 )
	return(Dinput(in.p1q,in.p5q,in.p3q,in.p7q,in.p2q,in.p6q,
		      in.m1q,in.m0q,in.m2q,in.m3q));
	//return(Dinput(-in.p1,in.p5,in.p3,in.m1q,in.m0q,in.m2q,in.m3q));
      if (fabs( in.p6q - in.m3q ) < deltam2)
	return(Dinput(in.p1q,in.p6q,in.p3q,in.p2q,in.p7q,in.p5q,
		 in.m1q,in.m0q,in.m3q,in.m2q));
	//return(Dinput(-in.p1,in.p6,-in.p3,in.m1q,in.m0q,in.m3q,in.m2q));
    }

  cout << "fixM2a: Something is wrong: tried to move zero"
       << " mass leg to pos. 3, \n "
       << " but there is no zero mass leg" << endl;
  in.printmaple();
  exit(1);
}

Dinput fixM2bc(const Dinput & in){

  if ( (in.m0q != 0.0) || (in.m1q != 0.0) ){
    cout << "fixM2bc: m0,m1 should be already zero, exit program";
    exit(1);
  } 

  if ((fabs( in.p2q - in.m3q ) < deltam2 ) 
      && (fabs( in.p2q - in.m3q ) < deltam2 ))
    return(in);

  if ((fabs( in.p5q - in.m3q ) < deltam2 ) 
      && (fabs( in.p7q - in.m3q ) < deltam2 ))
    return(Dinput(in.p1q,in.p7q,in.p3q,in.p5q,in.p6q,in.p2q,
		  in.m0q,in.m1q,in.m3q,in.m2q));
    //return(Dinput(in.p1,in.p7,-in.p3,in.m0q,in.m1q,in.m3q,in.m2q));

  if ((fabs( in.p2q - in.m3q ) < deltam2 ) 
      && ( fabs( in.p6q - in.m3q ) < deltam2 ) )
    return(Dinput(in.p1q,in.p6q,in.p3q,in.p2q,in.p7q,in.p5q,
		  in.m1q,in.m0q,in.m3q,in.m2q));
    //return(Dinput(-in.p1,in.p6,-in.p3,in.m1q,in.m0q,in.m3q,in.m2q));
 
  cout << "fixM2bc: Something is wrong: tried to move zero"
       << " mass leg to pos. 3, \n "
       << " but there is no zero mass leg" << endl;
  in.printmaple();
  exit(1);
}



/*********************************************************************
 *********************************************************************
 ***                                    ******************************
 *** Evaluation of the scalar integrals ******************************
 ***                                    ******************************
 *********************************************************************
 *********************************************************************
 *********************************************************************/

double beta(double r,double mq){
  return(sqrt(1.0-4.0*mq/r));
}

double xr(double beta){
  return((beta-1.0)/(beta+1.0));
}

complex<double> cclog1(complex<double> z){

  double zr = z.real();
  if ( zr > 0.0 ) 
    return complex<double>(log(zr),0.0);
  else 
    return complex<double>(log(-zr), z.imag() > 0.0 ? M_PI : -M_PI); 
}

complex<double> cclog2(double z, double ieps){

  if ( z > 0.0 ) 
    return complex<double>(log(z),0.0);
  else 
    return complex<double>(log(-z), ieps > 0.0 ? M_PI : -M_PI); 
}

complex<double> ccLi22(double x,double ieps) {
  if (x<1.0)
    return complex<double>(li2_(x),0.0);
  else
    return complex<double>(li2_(x),ieps > 0 ? M_PI * log(x) : - M_PI * log(x));
}

complex<double> ccdilog(double s1, double s2){
  /*
   * Evaluation of Li2(1-s1/s2) it is assumed that both s1 and s2 
   * have a small positive imaginary part.
   */
  if (s1*s2>0.0) 
    return complex<double>(li2_(1.0-s1/s2),0.0);
  else {
    double tmin = 1.0/(1.0-s1/s2);
    return complex<double>(li2_(1.0-s1/s2), s1<0.0 ? M_PI * log(tmin) : 
			   - M_PI * log(tmin) );
  }
}

complex<double> calLi2(double xr, double xi, double yr, double yi){

  if (xr*yr < 1.0){
    return(li2_(1.0-xr*yr) 
	   + log(1.0-xr*yr) * ( log(fabs(xr*yr)) - cclog2(xr,xi)
				- cclog2(yr,yi) ));
  } else {
    if (xi*yi<0.0) 
      return complex<double>(li2_(1.0-xr*yr),0.0);
    else
    cout << "calLi2: x*y > 1 and Im(x)*Im(y)>0 this case should never occur."
	 << " Exit...\n";
    exit(1);
  }
}



IntType I30(const Cinput & in){
  
  complex<double> res;

  /*
   * Check that all internal masses are zero
   */
  if ( in.m0q + in.m1q + in.m2q != 0.0 ) {
    cout << "I30: Error, non-zero internal mass...\n";
    in.printmaple();
    exit(1);
  }

  if (in.getExtMassCount() == 3){
    cout << "I30: All external legs are massive this case is not yet\n"
	 << "implemented. Exit...\n";
    in.printmaple();
    exit(1);
  }

  if (in.getExtMassCount() == 1){

    double s,
      s1 = in.p1q,
      s2 = in.p2q,
      s12 = in.p5q;

    if ( ( fabs(s1) > tiny ) && ( fabs(s2) < tiny ) && (fabs(s12) < tiny ) ) {
      s = s1;
    } else {
      if ( ( fabs(s2) > tiny ) && ( fabs(s1) < tiny ) 
	   && (fabs(s12) < tiny ) ) {
	s = s2;
      } else {
	if ( ( fabs(s12) > tiny ) && ( fabs(s1) < tiny ) 
	     && (fabs(s2) < tiny ) ) {
	  s = s12;
	} else {
	    cout << "I30: more than one massive leg...";
	    in.printmaple();
	    exit(1);	  
	}	
      }
    }

    /*
     * First integral in eq. A.2 [hep-ph/0211352]
     */
    res = 1.0/s 
      * ( parms.getDeltaIR2() 
	  - parms.getDeltaIR1() 
	  * cclog2(-s/parms.getInternalScaleSquared(),-1.0)
	  + 1.0/2.0 * pow(cclog2(-s/parms.getInternalScaleSquared(),-1.0),2)
	  - pi2/6.0 );

  } 

  if (in.getExtMassCount() == 2) {

    double s1,s2;

    if ( fabs(in.p1q) < tiny  ){
      s1 = in.p2q;
      s2 = in.p5q;
    } else {
      if ( fabs(in.p2q) < tiny ){
	s1 = in.p1q;
	s2 = in.p5q;
      } else {
	if ( fabs(in.p5q) < tiny ){
	  s1 = in.p1q;
	  s2 = in.p2q;
	} else {
	  cout << "I30: There is no massless leg, exit...\n";
	  in.printmaple();
	  exit(1);
	}
      }
    }

    if ( ( fabs(s1) < tiny ) || ( fabs(s2) < tiny ) ){
      in.printmaple();
      cout << "s1 = " << s1 << endl;
      cout << "s2 = " << s2 << endl;
    }

    res = 1.0 / ( s1-s2) *
      ( 
       - parms.getDeltaIR1() 
       * ( + cclog2(-s1/parms.getInternalScaleSquared(),-1.0) 
	   - cclog2(-s2/parms.getInternalScaleSquared(),-1.0))
       + 1.0/2.0 * ( + pow(cclog2(-s1/parms.getInternalScaleSquared(),-1.0),2) 
		     - pow(cclog2(-s2/parms.getInternalScaleSquared(),-1.0),2))
       );
    
  } 

  return mycast(res);
}

IntType I31(const Cinput & in) {

  double mq = in.m0q + in.m1q + in.m2q ;
  double s1,s2;

  if (fabs(in.p1q) < tiny) {
    s1 = in.p2q;
    s2 = in.p5q;
  } else 
    if (fabs(in.p2q) < tiny){
      s1 = in.p1q;
      s2 = in.p5q;
    } else 
      if (fabs(in.p5q) < tiny){
	s1 = in.p1q;
	s2 = in.p2q;
      } else {
	cout << "I31: Apparently this is not an I31 topology\n";
	in.printmaple();
	exit(1);
      }

  complex<double> res;  

  if ( ( fabs( s1 - mq ) > deltam2 ) && ( fabs( s2 - mq ) > deltam2 )){
    /*
     * Two off-shell legs,
     * use third formula from eq.A2 [hep-ph/0211352]  
     *
     */
    double s45 = s1;
    double t23 = s2;

    /*
     * If s45 == t23 we need a special formulae:
     */
    if ( fabs(s45-t23)/ ( fabs(s45) + fabs(t23) ) < 1.e-10 ) {
      res = -1.0 / ( t23 - mq )
	* ( parms.getDeltaIR1() + log(parms.getInternalScaleSquared()/mq) )
	+ cclog2(1.0-t23/mq,-1.0) * ( t23 + mq ) / ( t23 - mq ) / t23;
    } else {
      res = 1.0 / (s45-t23) 
	* ( + ( parms.getDeltaIR1() + log(parms.getInternalScaleSquared()/mq) )
	    * cclog2((t23-mq)/(s45-mq),s45-t23)
	    + ccLi22(s45/mq,1.0) - ccLi22(t23/mq,1.0)
	    + pow(cclog2(1.0-s45/mq,-1.0),2) 
	    - pow(cclog2(1.0-t23/mq,-1.0),2) );
    }
  } else {
    double t13;
    if (fabs(s1-mq) < deltam2)
      t13 = s2;
    else
      t13 = s1;
    /*
     * Use second formula from eq.A2 [hep-ph/0211352]  
     */

    res = 1.0/(t13-mq)
      * ( + 1.0/2.0 
	  * ( parms.getDeltaIR2() 
	      + parms.getDeltaIR1() 
	      * log(parms.getInternalScaleSquared()/mq)
	      + 1.0/2.0 * pow(log(parms.getInternalScaleSquared()/mq),2)
	      )
	  - ( parms.getDeltaIR1() + log(parms.getInternalScaleSquared()/mq) )
	  * cclog2(1.0-t13/mq,-1.0)
	  + ccLi22(t13/mq,1.0)
	  + pow(cclog2(1.0-t13/mq,-1),2));    
    
  }

  return(mycast(res));

}

IntType I32(double s34, double mtq){

  /*
   * Last integral in A.2 from hep-ph/0211352
   */
  
  double beta_s34 = beta(s34,mtq);//+ieps
  double xs34 = xr(beta_s34); //+ieps

  return(mycast(1.0/s34/beta_s34 
		* ( ( + parms.getDeltaIR1() 
		      + log(parms.getInternalScaleSquared()/mtq) ) * cclog2(xs34,1.0) 
		    - 2.0 * ccLi22(-xs34,-1.0)
		    - 2.0*cclog2(xs34,1.0)*cclog2(1.0+xs34,1.0)
		    + 1.0/2.0*pow(cclog2(xs34,1.0),2)
		    - pi2/6.0)
		)
	 );

}

IntType M0(const Dinput & in){
  /*
   *
   * All internal masses are zero and only one external leg has non-zero
   * p^2.
   *
   */

  // Check that this is indeed a M0-topology
  if (fabs(in.m0q)+fabs(in.m1q)+fabs(in.m2q)+fabs(in.m3q) != 0.0 ){
    cout << "M0: At least one of the masses seems to be a bit large\n"
	 << "to be considered as massless, exit programm.\n";
    in.printmaple();
    exit(1);
  }

  if ( fabs(in.p1q) + fabs(in.p2q) + fabs(in.p3q) > 1.E-8 ){
    cout << "M0: At least one of the external masses seems to be a bit large\n"
	 << "to be considered as massless, exit programm.\n";
    in.printmaple();
    exit(1);
  }

    
  double s12  = in.p5q;
  double s23  = in.p7q;
  double s123 = in.p6q;



  return(mycast( -2.0 /  s12 / s23 
		 * ( - parms.getDeltaIR2()
		     + parms.getDeltaIR1() * ( + cclog2(-s12/parms.getInternalScaleSquared(),-1.0)
				    + cclog2(-s23/parms.getInternalScaleSquared(),-1.0) 
				    - cclog2(-s123/parms.getInternalScaleSquared(),-1.0) ) 
		     - 0.5 * ( + pow(cclog2(-s12/parms.getInternalScaleSquared(),-1.0),2)  
			       + pow(cclog2(-s23/parms.getInternalScaleSquared(),-1.0),2)
			       - pow(cclog2(-s123/parms.getInternalScaleSquared(),-1.0),2) )  
		     + ccdilog(s123,s12) + ccdilog(s123,s23) 
		     + 1.0/2.0 * pow(cclog2(s12/s23,s23-s12),2) 
		     + pi2/3.0 
		     )
		 )
	 );

}

IntType M1a(const Dinput & in) {  
  /*
   * First integral in A.4 from [ehp-ph0211352]
   */
  double mtq = in.m3q;
  double shat = in.p5q;
  double t13 = in.p7q;
  double s45 = in.p6q;

  return(mycast( 
		1.0/shat/(t13-mtq) 
		* ( 3.0/2.0 * ( parms.getDeltaIR2() + parms.getDeltaIR1() * log(parms.getInternalScaleSquared()/mtq)
				+ 1.0/2.0 * pow(log(parms.getInternalScaleSquared()/mtq),2) 
				)
		    - (parms.getDeltaIR1() +  log(parms.getInternalScaleSquared()/mtq) )
		       * (2.0*cclog2(1.0-t13/mtq,-1.0) 
				+ cclog2(-shat/mtq,-1.0) 
				- cclog2(1.0-s45/mtq,-1.0))
		    - 2.0*ccdilog(s45-mtq,t13-mtq)
		    + 2.0 * cclog2(-shat/mtq,-1.0)
		    * cclog2(1.0-t13/mtq,-1.0)
		    - pow(cclog2(1.0-s45/mtq,-1.0),2)
		    - 2.0*pi2/3.0 )
		) 
	 );
}

IntType M1b(const Dinput & in) {  

  double mtq = in.m3q;
  double s = in.p5q;
  double t = in.p7q;
  double m_3q = in.p2q;

  return(mycast(1.0/s/(t-mtq) *
		( + 1.0/2.0 * ( parms.getDeltaIR2()
				+ parms.getDeltaIR1() * log(parms.getInternalScaleSquared()/mtq)
				+ 1.0/2.0 * pow(log(parms.getInternalScaleSquared()/mtq),2) 
				)
		  - ( parms.getDeltaIR1() + log(parms.getInternalScaleSquared()/mtq) )
		  * ( + cclog2(1.0-t/mtq,-1.0)
		      + cclog2(-s/mtq,-1.0)
		      - cclog2(-m_3q/mtq,-1.0))
		  + 1.0/2.0*pow(cclog2(1.0-t/mtq,-1.0),2)
		  + 1.0/2.0*pow(cclog2(-s/mtq,-1.0),2)
		  - 1.0/2.0*pow(cclog2(-m_3q/mtq,-1.0),2)
		  + 2.0 * ccdilog(s,m_3q)
		  + 1.0/2.0*pow(cclog2(s/m_3q,m_3q-s)
				+ 2.0*cclog2(1.0-t/mtq,-1.0),2)
		  - cclog2(-m_3q/mtq,-1.0)*cclog2(s/m_3q,m_3q-s)
		  - 3.0/2.0*pow(cclog2(1.0-t/mtq,-1.0),2)
		  + pi2/6.0)
		)
	 );
}

IntType M2a(const Dinput & in) {  

  /*
   * Second integral in A.4 from [ehp-ph0211352]
   */
  double t13 = in.p5q;
  double t23 = in.p7q;
  double s45 = in.p6q;
  
  double mtq = in.m2q;

  return(mycast(1.0/(t13-mtq)/(t23-mtq)
		* ( 1.0/2.0 * ( parms.getDeltaIR2() 
				+ parms.getDeltaIR1() * log(parms.getInternalScaleSquared()/mtq)
				+ 1.0/2.0 * pow(log(parms.getInternalScaleSquared()/mtq),2) 
				)
		    - ( parms.getDeltaIR1() + log(parms.getInternalScaleSquared()/mtq) ) 
		    * ( + cclog2(1.0 - t13/mtq,-1.0)
			+ cclog2(1.0 - t23/mtq,-1.0)
			- cclog2(1.0 - s45/mtq,-1.0) )
		    - 2.0*ccdilog(s45-mtq,t23-mtq)
		    - 2.0*ccdilog(s45-mtq,t13-mtq)
		    + 2.0*cclog2(1.0 - t13/mtq,-1.0) * cclog2(1.0 - t23/mtq,-1.0)
		    - pow(cclog2(1.0 - s45/mtq,-1.0),2)
		    - pi2/6.0 )));
  

}

IntType M2b(const Dinput & in) {  

  /*
   * Third integral in A.4 from [ehp-ph0211352]
   */
  double t13 = in.p5q;
  double t14 = in.p7q;
  double t25 = in.p3q;
  double mtq = in.m2q;
  double beta_t25 = beta(t25,mtq); // + i eps
  double xt25 = xr(beta_t25); // + ieps

  return(mycast( 
		1.0/(t13-mtq)/(t14-mtq)
		* ( + ( parms.getDeltaIR2()	
			+ parms.getDeltaIR1() * log(parms.getInternalScaleSquared()/mtq)
			+ 1.0/2.0 * pow(log(parms.getInternalScaleSquared()/mtq),2) 
			)
		    - (parms.getDeltaIR1() + log(parms.getInternalScaleSquared()/mtq) )
		    *(+ cclog2(1.0-t14/mtq,-1.0)
		      + cclog2(1.0-t13/mtq,-1.0)
		      )
		    - pow(cclog2(xt25,1.0),2)
		    + 2.0*cclog2(1.0-t13/mtq,-1.0)
		    * cclog2(1.0-t14/mtq,-1.0)
		    - 2.0*pi2/3.0 ) 
		)
	 );

}


IntType M3a(const Dinput & in) {  

  /***
   *** Check that the input is correct:
   ***/
  if ( (fabs(in.p1q-in.p6q)>1.0) || (fabs(in.p1q-in.m1q) > deltam2 ) 
       || (in.m0q != 0.0) ){
    cout << "M3a: Input has not the appropriate format\n";
    in.printmaple();
    exit(1);
  }


  double t13 = in.p5q,
    s34 = in.p7q,
    t25 = in.p2q,
    mtq = in.m1q;

  double beta_s34 = beta(s34,mtq); // + i eps
  double beta_t25 = beta(t25,mtq); // + i eps
  double xs34 = xr(beta_s34); //+ieps
  double xt25 = xr(beta_t25); //+ieps
  
  return(mycast(1.0/(t13-mtq)/s34/beta_s34
		* ( ( parms.getDeltaIR1() + log(parms.getInternalScaleSquared()/mtq) )
		    * cclog2(xs34,1.0) 
		    - 2.0 * ( 
			     + calLi2(xs34,1.0,xt25,1.0) 
			     + calLi2(xs34,1.0,1.0/xt25,-1.0)
			     )
		    - ccLi22(xs34*xs34,xs34) 
		    - 2.0*cclog2(xs34,1.0)*cclog2(1.0-xs34*xs34,-xs34)
		    - 2.0*cclog2(xs34,1.0)*cclog2(1.0-t13/mtq,-1.0) 
		    - pow(cclog2(xt25,1.0),2) 
		    + pi2/6.0)
		)
	 ); 

}

IntType ScalarInt(const double mq){

  //  return mycast(complex<double>(parms.getDeltaUV1(),0.0));
  
  if (!initialized) initIntegralLibrary();
  
  if ( mq == 0.0 ) {
    return mycast(complex<double>(mq*(parms.getDeltaUV1()-parms.getDeltaIR1()),0.0));
    //return mycast(complex<double>((parms.getDeltaUV1()-parms.getDeltaIR1()),0.0)); // 06.06.06
  } else {
    //        return mycast(complex<double>(1.e30,0.0));
    if (OnlyUVSing) return mycast(complex<double>(mq,0.0));
    
    /*
     * Use FF lib to calculate the finite integral:
     */
    complex<double> cint;
    int ierr=0;    
    ffxa0_(cint,parms.getDeltaUV1(),parms.getInternalScaleSquared(),
	   mq,ierr);
    return mycast(cint);
  }
}

#define returnA(__XX__) {double tmp = __XX__;\
cout <<" A ( "<< mq << " ) = " << tmp << endl;\
return (tmp);}

#define returnB(__XX__) {double tmp = __XX__;\
cout <<" B = " << tmp << endl;\
return (tmp);}

#define returnC(__XX__) {double tmp = __XX__;\
cout <<" C = " << tmp << endl;\
return (tmp);}

#define returndB(__XX__) {double tmp = __XX__;\
cout <<" diffB = " << tmp << endl;\
return (tmp);}

IntType A(double mq){
  return(ScalarInt(mq));
}

IntType B(double pq, double m0q, double m1q){
  /*
  complex<double> ffresult = ScalarInt(pq, m0q, m1q);
  complex<double> denner = b0check_(parms.getDeltaUV1(), 
		   parms.getRenormalizationScaleSquared(),
				    pq,m0q,m1q);
  
  cout.precision(15);
  cout << mycast(denner/ffresult) << endl;
  */
  return(ScalarInt(pq, m0q, m1q));
}

IntType C(double p1q, double p2q, double p5q, 
	  double m0q,double m1q, double m2q){
  return( ScalarInt( Cinput(p1q,p2q,p5q,m0q,m1q,m2q) ) );
}

IntType D(double p1q, double p2q, double p3q, double p4q, 
	  double p5q, double p7q, 
	  double m0q,double m1q, double m2q, double m3q){
  return( ScalarInt( Dinput(p1q,p2q,p3q,p4q,p5q,p7q,m0q,m1q,m2q,m3q) ) );
}

IntType diffB0(double pq, double m0, double m1){
  complex<double> r,cint;

  if ( m0*m1 == 0.0 ) {
    cout << "wrong arguments in RediffB0, m0, m1 must be non-zero\n";
    exit(1);
  }
  double m0q = m0*m0;
  double m1q = m1*m1;
  r = ( m0q+m1q-pq
       + sqrt( complex<double>( (-pq+m1q-2.0*m0*m1+m0q) * 
				(-pq+m1q+2.0*m0*m1+m0q) ) ) ) /m0/m1/2.0;

  if ( r == complex<double>(1.0,0.0) ) 
    cint = - (m0q-m1q)/pq/pq*log(m1/m0) - 1.0/pq*2.0;
  else   
    cint = - (m0q-m1q)/pq/pq*log(m1/m0)
      + m0*m1/pq/pq*(1.0/r-r)*log(r)
      - 1.0/pq*(1.0+(r*r+1.0)/(r*r-1.0)*log(r));
  return(mycast(cint));
}

IntType ScalarInt(double p1q, double m0q, double m1q) {

  //  return mycast(complex<double>(0*parms.getDeltaUV1(),0.0));
  if (!initialized) initIntegralLibrary();

  if ( (fabs(p1q) < 1E-7) && (m0q < 1E-8) && (m1q < 1E-8) ) {
    //cout << "ScalarInt(B): How could we end up here?\n";
    return mycast(complex<double>(parms.getDeltaUV1()-parms.getDeltaIR1(),0.0));
  } 

  if (  (fabs(p1q) < 1E-7) && (m0q==m1q) ) {
    return(ScalarInt(m0q)/m0q-1.0);
  }
  if (OnlyUVSing) return mycast(complex<double>(0.0,0.0));

  /*
   * Use FF lib to calculate the finite integral:
   */
  complex<double> cint;
  int ierr=0;    
  ffxb0_(cint, parms.getDeltaUV1(), parms.getInternalScaleSquared(), 
	 p1q , m0q, m1q, ierr);
  return mycast(cint);
}

IntType ScalarInt(const Binput & in){ 

  return(ScalarInt(in.p1q,in.m0q,in.m1q));
}

IntType ScalarInt(const Cinput & in){

  //  return mycast(complex<double>(1e20,0.0));

  if (!initialized) initIntegralLibrary();
  if (OnlyUVSing) return mycast(complex<double>(0.0,0.0));

  switch (in.getMassCount()){
  case 0:
    {
      /*
      * All internal masses are zero:
      */ 
      return(I30(in));
    }
    break;
  case 1:
    { 
      /*
       * Two internal masses are zero, one is non-zero
       */
      // check if there is one massless external leg:
      if ( ( fabs(in.p1q) < tiny ) || ( fabs(in.p2q) < tiny )
	   || ( fabs(in.p5q) < tiny ) ) {
	return(I31(in));
      }

      /*
       * Let FF do the job...
       */
    }
    break;
  case 2:
    {
      /*
      * Two internal masses --> one IR divergent triangle
      */
      double mq = (in.m0q+in.m1q+in.m2q) / 2.0;

      if ( ( fabs( in.p1q - mq ) < deltam2 ) &&
	   ( fabs( in.p2q - mq ) < deltam2 ) )
	return( I32(in.p5q,mq) );
      
      if ( ( fabs( in.p1q - mq ) < deltam2 ) &&
	   ( fabs( in.p5q - mq ) < deltam2 ) )
	return( I32(in.p2q,mq) );

      if ( ( fabs( in.p2q - mq ) < deltam2 ) &&
	   ( fabs( in.p5q - mq ) < deltam2 ) )
	return( I32(in.p1q,mq) );
      /*
       * Not divergent, use FF...
       */
    }
    break;
  case 3:
    /*
     * All internal masses non-zero, let FF do the job...
     */
    break;
  }

  /*
   * Use FF lib to calculate the finite integral:
   */

  complex <double> cint;
  double xxpi[6];
  int ierr=0;
  xxpi[0] = in.m0q;
  xxpi[1] = in.m1q;
  xxpi[2] = in.m2q;
  xxpi[3] = in.p1q; 
  xxpi[4] = in.p2q; 
  xxpi[5] = in.p5q; 
  ffxc0_(cint,xxpi,ierr);
  return mycast( cint);
}


IntType ScalarInt(const Dinput & in){

  //  return mycast(complex<double>(0e20,0.0));

  if (!initialized) initIntegralLibrary();
  if (OnlyUVSing) return(mycast(complex<double>(0.0,0.0)));

  /*
   * Figure out the exact topology:
   */
  switch (in.getMassCount()){
  case 0:
    { /*
       * All internal masses are zero
       */
      Dinput tmp=fixM0(in);
      return(M0(tmp));
    }
    break;
  case 1:
    {
      /*
       * Three internal masses are zero, one is non zero.
       *
       * Figure out how many external legs are massless...
       */
      int ext = 6-in.getExtMassCount();
      switch (ext) {
      case 2: {
	/*
	 * Two external legs are massive
	 */
	Dinput tmp=fixM1a2(fixM1a1(fixDm123(in)));
	return(M1a(tmp));
      } break;
      case 1: {
	/*
	 * Three external legs are massive
	 */
	Dinput tmp=fixM1b2(fixM1b1(fixDm123(in)));
	return(M1b(tmp));
      } break;
      }
    }
    break;
  case 2:
    {
      /*
       * Two internal masses are zero,
       *
       * figure out how many external legs are massless
       */
      int ext = 6-in.getExtMassCount();
      switch (ext) {
      case 2: {
	Dinput tmp=fixM2a(fixDm1(fixDm0(in)));
	return(M2a(tmp));
      } break;
      case 1: {
	Dinput tmp=fixM2bc(fixDm1(fixDm0(in)));
	if ( fabs(tmp.p1q) < tiny ){
	  return(M2b(tmp));
	}
	if ( fabs(tmp.p3q) < tiny ){
	  // M2c do nothing, use FF lib to evaluate this integral
	} else {
	  cout << "Something went wrong exit...\n";
	  exit(1);
	}
      } break;
      }
    }
    break;
  case 3:
    { 
      int ext = 6-in.getExtMassCount();
      switch (ext) {
      case 1:{
	Dinput tmp = fixM3a2(fixM3a1(fixDm0(in))); 
	return(M3a(tmp));
      } break;
      case 2: {
	// Do nothing, this integral is evaluated using the FF lib.
	// Dinput tmp=fixM3b(fixDm0(in));
	// return(M3b(tmp));
      } break;
      }
    }
    break;
    case 4:
      /*
       * All internal masses non-zero, this case is handled by FF.
       */
    break;
  }
  /*
   * Use FF lib to calculate the finite integral:
   */
  complex <double> cint;
  double xxpi[13];
  int ierr = 0;
  xxpi[0] = in.m0q;
  xxpi[1] = in.m1q;
  xxpi[2] = in.m2q;
  xxpi[3] = in.m3q;
  xxpi[4] = in.p1q; 
  xxpi[5] = in.p2q; 
  xxpi[6] = in.p3q; 
  xxpi[7] = in.p6q; 
  xxpi[8] = in.p5q;
  xxpi[9] = in.p7q;
  xxpi[10] = 0.0;
  xxpi[11] = 0.0;
  xxpi[12] = 0.0;  
  ffxd0_(cint,xxpi,ierr);
  //if (ierr > 10 )
  //in.printmaple();
  return(mycast(cint));
}


IntType ScalarInt(const vector<FourMomentum>& q, const vector<double>& m){

  if ((q.size() == 2) && (m.size() == 2)){
    return(ScalarInt(Binput(dotp(q[1]-q[0],q[1]-q[0]),m[0]*m[0],m[1]*m[1])));
  } 
  if ((q.size() == 3) && (m.size() == 3)){
    return(ScalarInt(Cinput(dotp(q[1]-q[0],q[1]-q[0]),
			    dotp(q[2]-q[1],q[2]-q[1]),
			    dotp(q[2]-q[0],q[2]-q[0]),
		  m[0]*m[0],m[1]*m[1],m[2]*m[2])));
  } 
  if ((q.size() == 4) && (m.size() == 4)){
    //    return(ScalarInt(Dinput(q[1]-q[0],q[2]-q[1]-q[0],q[3]-q[2]-q[1]-q[0],
    //   data help(dbox(q[0],q[1],q[2],q[3],m[0],m[1],m[2],m[3]));
    //cout.precision(15);
    //cout <<"Determinant: " << help.detS << endl;
    //    return(ScalarInt(Dinput(q[1]-q[0],q[2]-q[1],q[3]-q[2],
    //		  m[0]*m[0],m[1]*m[1],m[2]*m[2],m[3]*m[3])));
    return(ScalarInt(Dinput(dotp(q[1]-q[0],q[1]-q[0]),
			    dotp(q[2]-q[1],q[2]-q[1]),
			    dotp(q[3]-q[2],q[3]-q[2]),
			    dotp(q[3]-q[0],q[3]-q[0]),
			    dotp(q[2]-q[0],q[2]-q[0]),
			    dotp(q[3]-q[1],q[3]-q[1]),
			    m[0]*m[0],m[1]*m[1],m[2]*m[2],m[3]*m[3])));
  }
  cout << " Cannot handle scalar integrals with " << q.size() 
       << "momenta and " << m.size() << " masses." << endl;
  exit(1);
  return(mycast(complex<double>(0.0,0.0)));
}

