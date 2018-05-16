// $Modified: Fri Jan 29 10:19:53 2010 by uwer $
#define NDEBUG

#include "StandardModelParameters.h"
#include "GGReduction.h"
#include "ScalarInt.h"
#include "zlib.h"
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <cassert>
#include <map>


#undef NEWREC
//#define NEWREC


#ifndef ONLYREAL
#define INT2DOUBLE(_X_) static_cast<double>(_X_)
#else
#define INT2DOUBLE(_X_) (_X_)
#endif

using namespace std;

static StandardModelParameters& parms = StandardModelParameters::instance();

const double expflag = 1.0;

static double BCUT = 5.e-13;

int GGException::count = 0;

void GGException::print(){

  if ( ( count % 10000 ) == 1 ){
    cout << "# " << what() << endl;
    if ( count != 1 ) 
      cout << "# Last message repeated 10000 times (total = " 
	   << count << " )" << endl;
  }
}

/*
 * Values of the Gammafunction for positive integer values of the arguments
 * 
 * gam[i] gives Gamma(i)
 *
 */
static const double gam[21] = {0.0,1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,
			 40320.0,362880.0,3628800.0,39916800.0,479001600.0,
			 6227020800.0,87178291200.0,1307674368000.0,
			 20922789888000.0,355687428096000.0,
			 6402373705728000.0,1.21645100408832000e+17}; 

/*
 * S1(n) = sum(1/i,i=1..n)
 * 
 */
static const double S1[31] = {0.0, 
			      1.0, 3.0/2.0, 11.0/6.0, 25.0/12.0,
			      137.0/60.0, 49.0/20.0, 363.0/140.0, 761.0/280.0, 
			      7129.0/2520.0, 7381.0/2520.0,83711.0/27720.0,
			      86021.0/27720.0,1145993.0/360360.0,
			      1171733.0/360360.0,1195757.0/360360.0,
			      2436559.0/720720.0,42142223.0/12252240.0,
			      14274301.0/4084080.0,275295799.0/77597520.0,
			      55835135.0/15519504.0,18858053.0/5173168.0,
			      19093197.0/5173168.0,444316699.0/118982864.0,
			      1347822955.0/356948592.0,
			      34052522467.0/8923714800.0,
			      34395742267.0/8923714800.0,
			      312536252003.0/80313433200.0,
			      315404588903.0/80313433200.0,
			      9227046511387.0/2329089562800.0,
			      9304682830147.0/2329089562800.0};

void setBCUT(double value){
  BCUT = value;
  cout << "# Set BCUT = " << BCUT << endl;
}

double setbcut(double value){
  double tmp = BCUT;
  BCUT = value;
  return(tmp);
}

/*
 * The following functions produces a hash value to identify integrals
 * with raised dimensions and raised propagators.
 *
 * At a certain point I should change the hash function:
 * - for ni > 10  we might get wrong results
 * - using a the basis a power of 2 might be much faster...   
 */
inline unsigned int ihash(unsigned int d, unsigned int n1, unsigned int n2, 
			  unsigned int n3, unsigned int n4, unsigned int n5,
			  unsigned int n6){
  //  return(d*100000000+n1*10000000+n2*1000000+n3*100000+n4*10000+n5*1000);
  return( (d << 24) + (n1 << 20 ) + ( n2 << 16 ) + ( n3 << 12 ) + (n4 << 8 ) + (n5 << 4 ) + n6);
}

inline unsigned int ihash(unsigned int d, unsigned int n1, unsigned int n2, 
		   unsigned int n3, unsigned int n4, unsigned int n5){
  //  return(d*100000000+n1*10000000+n2*1000000+n3*100000+n4*10000+n5*1000);
  return( (d << 24) + (n1 << 20 ) + ( n2 << 16 ) + ( n3 << 12 ) + (n4 << 8 ) + (n5 << 4 ));
}

inline unsigned ihash(unsigned int d, unsigned int n1, unsigned int n2, 
		      unsigned int n3, unsigned int n4){
  //return(d*100000000+n1*10000000+n2*1000000+n3*100000+n4*10000);
  return( (d << 24) + (n1 << 20 ) + ( n2 << 16 ) + ( n3 << 12 ) + (n4 << 8 ) );
}

inline unsigned ihash(unsigned int d, unsigned int n1, unsigned int n2, 
		      unsigned int n3){
  //  return(d*100000000+n1*10000000+n2*1000000+n3*100000);
  return( (d << 24) + (n1 << 20 ) + ( n2 << 16 ) + ( n3 << 12 ) );
}



/*
 *
 * This is the place where we store all the information about topologies
 * we have already seen, it is a static member of Topology and as such
 * needs to be initialized outside.
 *
 */
map<uLong,data*> data::topologies;


#ifdef WITHCACHE
/*
 * The following return commands are used to register a result
 * in the integral cache of the corresponding topology
 * before returning the value.
 */

#define RETURN6(VALUE) {IntType tmpxxx = (VALUE);\
(topo->cache)[ihash(d,n1,n2,n3,n4,n5,n6)] = tmpxxx;\
return(tmpxxx);\
}

#define RETURN5(VALUE) {IntType tmpxxx = (VALUE);\
(topo->cache)[ihash(d,n1,n2,n3,n4,n5)] = tmpxxx;\
return(tmpxxx);\
}

#define RETURN4(VALUE) {IntType tmpxxx = (VALUE);\
(topo->cache)[ihash(d,n1,n2,n3,n4)] = tmpxxx;\
return(tmpxxx);\
}

#define RETURN3(VALUE) {IntType tmpxxx = (VALUE);\
(topo->cache)[ihash(d,n1,n2,n3)] = tmpxxx;\
return(tmpxxx);\
}

#else 

/*
 * Don't use any cache, every integral is always evaluated by the
 * full recursion. Very slow, mainly for testing.
 */


#define RETURN6(x__) return(x__);
#define RETURN5(x__) return(x__);
#define RETURN4(x__) return(x__);
#define RETURN3(x__) return(x__);

#endif // ifdef WITHCACHE



#define PRINT(X) cout << #X <<" = " << X << endl


static inline int imax(const int n1, const int n2, const int n3, const int n4, 
	 const int n5, const int n6){
  int i=1,t=n1;
  i = (t > n2) ? i : (t=n2,2) ; 
  i = (t > n3) ? i : (t=n3,3) ;
  i = (t > n4) ? i : (t=n4,4) ;
  i = (t > n5) ? i : (t=n5,5) ;
  i = (t > n6) ? i : (t=n6,6) ;
  return(i);
} 

/******************************************************************
 *
 * Reduction of the 5-point topologies
 *
 *
 ******************************************************************/

IntType XInt( int d,  int n1,  int n2,  int n3, 
	      int n4,  int n5, Topology* topo, int order){

  cout << "# new5(" << d << "," << n1 << "," << n2 << ","<< n3 << "," 
       << n4 << "," << n5 <<")" << endl;

  if ( order > 1 )
    return ( 0.);

  const  int n = 5;
  const double B = topo->getB();
  const  int sigma = n1 + n2 + n3 + n4 + n5;
  double b[n+1];

  cout << "# order = " << order << endl; 
  for ( int i=1; i<= n; i++) {
    b[i] = 0.;
    for ( int j=1; j<= n; j++) {
      b[i] += topo->getSinv(i,j);
    }
  }

  return(//( d - 1.0 - sigma ) * B * Int(d+2,n1,n2,n3,n4,n5,topo,order+1)
	 + b[1] * Int(d,n1-1,n2,n3,n4,n5,topo)
	 + b[2] * Int(d,n1,n2-1,n3,n4,n5,topo)
	 + b[3] * Int(d,n1,n2,n3-1,n4,n5,topo)
	 + b[4] * Int(d,n1,n2,n3,n4-1,n5,topo)
	 + b[5] * Int(d,n1,n2,n3,n4,n5-1,topo) );

}

IntType XInt( int d,  int n1,  int n2,  int n3,
	     Topology* topo, int order){

  cout << "# new5(" << d << "," << n1 << "," << n2 << ","<< n3 
       <<")" << endl;

  if ( order > 1 )
    return ( 0.);

  const  int n = 3;
  const double B = topo->getB();
  const  int sigma = n1 + n2 + n3;
  double b[n+1];

  cout << "# order = " << order << endl; 
  for ( int i=1; i<= n; i++) {
    b[i] = 0.;
    for ( int j=1; j<= n; j++) {
      b[i] += topo->getSinv(i,j);
    }
  }

  return(//( d - 1.0 - sigma ) * B * Int(d+2,n1,n2,n3,topo,order+1)
	 + b[1] * Int(d,n1-1,n2,n3,topo)
	 + b[2] * Int(d,n1,n2-1,n3,topo)
	 + b[3] * Int(d,n1,n2,n3-1,topo)
	 );

}

/*********************************************************************
 * SIX-POINT REDUCTION
 *********************************************************************/
#define VERYNEW
#ifdef VERYNEW
IntType Int( int d,  int n1,  int n2,  int n3, 
	     int n4,  int n5, int n6, Topology* topo){


  if ( topo == 0 ) {
    cout << "# GGReduction: Null topology encountered" << endl;
    /***
     *** For testing purposes: return zero for undefined topologies.
     ***/ 
    //    return(0.);
    exit(1);
  }

#ifdef GGDEBUG
  cout << "# Hexagon(" << d << "," << n1 << "," << n2 << ","<< n3 << "," 
       << n4 << "," << n5 << n6 <<")" << endl;
#endif

#ifdef WITHCACHE

  /*
   * Check in the integral cache if we know this integral already:
   */

  map<unsigned int,IntType>::iterator integral = 
    (topo->cache).find(ihash(d,n1,n2,n3,n4,n5,n6));
  
  if ( integral != (topo->cache).end() ) {
    return( (*integral).second );
  }

#endif


  const  int n = 6;
  const  int sigma = n1 + n2 + n3 + n4 + n5 + n6;

  double B, b[n+1], Sinv[n+1][n+1];

  B = topo->getB();


  
  
  if (n1==0) RETURN6(Int(d,n2,n3,n4,n5,n6,topo->subtopo(1)));
  if (n2==0) RETURN6(Int(d,n1,n3,n4,n5,n6,topo->subtopo(2)));
  if (n3==0) RETURN6(Int(d,n1,n2,n4,n5,n6,topo->subtopo(3)));
  if (n4==0) RETURN6(Int(d,n1,n2,n3,n5,n6,topo->subtopo(4)));
  if (n5==0) RETURN6(Int(d,n1,n2,n3,n4,n6,topo->subtopo(5)));
  if (n6==0) RETURN6(Int(d,n1,n2,n3,n4,n5,topo->subtopo(6)));


  for (int i=1; i<= n; i++) {
    b[i] = 0.;
    for (int j=1; j<= n; j++) {
      Sinv[i][j] = topo->getSinv(i,j);
      b[i] += Sinv[i][j];
    }
  }


  if ( sigma == n || (d/2-sigma + 1 < 0) )
    RETURN6( + b[1] * Int(d,n1-1,n2,n3,n4,n5,n6,topo)
	     + b[2] * Int(d,n1,n2-1,n3,n4,n5,n6,topo)
	     + b[3] * Int(d,n1,n2,n3-1,n4,n5,n6,topo)
	     + b[4] * Int(d,n1,n2,n3,n4-1,n5,n6,topo)
	     + b[5] * Int(d,n1,n2,n3,n4,n5-1,n6,topo)
	     + b[6] * Int(d,n1,n2,n3,n4,n5,n6-1,topo)
	     );

  if ( 0 == d/2 - sigma + 1 ) {

    int k = imax(n1,n2,n3,n4,n5,n6);

    if ( 1 == k)
      RETURN6( -( INT2DOUBLE(n1) * ( + b[2] * Int(d,n1+1,n2-1,n3,n4,n5,n6,topo)
			 + b[3] * Int(d,n1+1,n2,n3-1,n4,n5,n6,topo)
			 + b[4] * Int(d,n1+1,n2,n3,n4-1,n5,n6,topo)
			 + b[5] * Int(d,n1+1,n2,n3,n4,n5-1,n6,topo)
			 + b[6] * Int(d,n1+1,n2,n3,n4,n5,n6-1,topo))
		  + Sinv[k][1] * Int(d-2,n1-1,n2,n3,n4,n5,n6,topo)
		  + Sinv[k][2] * Int(d-2,n1,n2-1,n3,n4,n5,n6,topo)
		  + Sinv[k][3] * Int(d-2,n1,n2,n3-1,n4,n5,n6,topo)
		  + Sinv[k][4] * Int(d-2,n1,n2,n3,n4-1,n5,n6,topo)
		  + Sinv[k][5] * Int(d-2,n1,n2,n3,n4,n5-1,n6,topo)
		  + Sinv[k][6] * Int(d-2,n1,n2,n3,n4,n5,n6-1,topo) ) 
	       / ((d-sigma+n1-1) * b[1]) );

    if ( 2 == k)
      RETURN6( -( INT2DOUBLE(n2) * ( + b[1] * Int(d,n1-1,n2+1,n3,n4,n5,n6,topo)
			 + b[3] * Int(d,n1,n2+1,n3-1,n4,n5,n6,topo)
			 + b[4] * Int(d,n1,n2+1,n3,n4-1,n5,n6,topo)
			 + b[5] * Int(d,n1,n2+1,n3,n4,n5-1,n6,topo)
			 + b[6] * Int(d,n1,n2+1,n3,n4,n5,n6-1,topo))
		  + Sinv[k][1] * Int(d-2,n1-1,n2,n3,n4,n5,n6,topo)
		  + Sinv[k][2] * Int(d-2,n1,n2-1,n3,n4,n5,n6,topo)
		  + Sinv[k][3] * Int(d-2,n1,n2,n3-1,n4,n5,n6,topo)
		  + Sinv[k][4] * Int(d-2,n1,n2,n3,n4-1,n5,n6,topo)
		  + Sinv[k][5] * Int(d-2,n1,n2,n3,n4,n5-1,n6,topo)
		  + Sinv[k][6] * Int(d-2,n1,n2,n3,n4,n5,n6-1,topo) ) 
	       / ((d-sigma+n2-1) * b[2]) );

      if ( 3 == k)
	RETURN6( -( INT2DOUBLE(n3) * ( + b[1] * Int(d,n1-1,n2,n3+1,n4,n5,n6,topo)
			 + b[2] * Int(d,n1,n2-1,n3+1,n4,n5,n6,topo)
			 + b[4] * Int(d,n1,n2,n3+1,n4-1,n5,n6,topo)
			 + b[5] * Int(d,n1,n2,n3+1,n4,n5-1,n6,topo)
			 + b[6] * Int(d,n1,n2,n3+1,n4,n5,n6-1,topo))
		  + Sinv[k][1] * Int(d-2,n1-1,n2,n3,n4,n5,n6,topo)
		  + Sinv[k][2] * Int(d-2,n1,n2-1,n3,n4,n5,n6,topo)
		  + Sinv[k][3] * Int(d-2,n1,n2,n3-1,n4,n5,n6,topo)
		  + Sinv[k][4] * Int(d-2,n1,n2,n3,n4-1,n5,n6,topo)
		  + Sinv[k][5] * Int(d-2,n1,n2,n3,n4,n5-1,n6,topo)
		  + Sinv[k][6] * Int(d-2,n1,n2,n3,n4,n5,n6-1,topo) ) 
		 / ((d-sigma+n3-1) * b[3] ) );

      if ( 4 == k)
	RETURN6( -( INT2DOUBLE(n4) * ( + b[1] * Int(d,n1-1,n2,n3,n4+1,n5,n6,topo)
			 + b[2] * Int(d,n1,n2-1,n3,n4+1,n5,n6,topo)
			 + b[3] * Int(d,n1,n2,n3-1,n4+1,n5,n6,topo)
			 + b[5] * Int(d,n1,n2,n3,n4+1,n5-1,n6,topo)
			 + b[6] * Int(d,n1,n2,n3,n4+1,n5,n6-1,topo))
		  + Sinv[k][1] * Int(d-2,n1-1,n2,n3,n4,n5,n6,topo)
		  + Sinv[k][2] * Int(d-2,n1,n2-1,n3,n4,n5,n6,topo)
		  + Sinv[k][3] * Int(d-2,n1,n2,n3-1,n4,n5,n6,topo)
		  + Sinv[k][4] * Int(d-2,n1,n2,n3,n4-1,n5,n6,topo)
		  + Sinv[k][5] * Int(d-2,n1,n2,n3,n4,n5-1,n6,topo)
		  + Sinv[k][6] * Int(d-2,n1,n2,n3,n4,n5,n6-1,topo) ) 
		 / ((d-sigma+n4-1) * b[4]) );

      if ( 5 == k)
	RETURN6( -( INT2DOUBLE(n5) * ( + b[1] * Int(d,n1-1,n2,n3,n4,n5+1,n6,topo)
			 + b[2] * Int(d,n1,n2-1,n3,n4,n5+1,n6,topo)
			 + b[3] * Int(d,n1,n2,n3-1,n4,n5+1,n6,topo)
			 + b[4] * Int(d,n1,n2,n3,n4-1,n5+1,n6,topo)
			 + b[6] * Int(d,n1,n2,n3,n4,n5+1,n6-1,topo))
		  + Sinv[k][1] * Int(d-2,n1-1,n2,n3,n4,n5,n6,topo)
		  + Sinv[k][2] * Int(d-2,n1,n2-1,n3,n4,n5,n6,topo)
		  + Sinv[k][3] * Int(d-2,n1,n2,n3-1,n4,n5,n6,topo)
		  + Sinv[k][4] * Int(d-2,n1,n2,n3,n4-1,n5,n6,topo)
		  + Sinv[k][5] * Int(d-2,n1,n2,n3,n4,n5-1,n6,topo)
		  + Sinv[k][6] * Int(d-2,n1,n2,n3,n4,n5,n6-1,topo) ) 
		 / ((d-sigma+n5-1) * b[5]) );

      if ( 6 == k)
	RETURN6( -( INT2DOUBLE(n6) * ( + b[1] * Int(d,n1-1,n2,n3,n4,n5,n6+1,topo)
			 + b[2] * Int(d,n1,n2-1,n3,n4,n5,n6+1,topo)
			 + b[3] * Int(d,n1,n2,n3-1,n4,n5,n6+1,topo)
			 + b[4] * Int(d,n1,n2,n3,n4-1,n5,n6+1,topo)
			 + b[5] * Int(d,n1,n2,n3,n4,n5-1,n6+1,topo))
		  + Sinv[k][1] * Int(d-2,n1-1,n2,n3,n4,n5,n6,topo)
		  + Sinv[k][2] * Int(d-2,n1,n2-1,n3,n4,n5,n6,topo)
		  + Sinv[k][3] * Int(d-2,n1,n2,n3-1,n4,n5,n6,topo)
		  + Sinv[k][4] * Int(d-2,n1,n2,n3,n4-1,n5,n6,topo)
		  + Sinv[k][5] * Int(d-2,n1,n2,n3,n4,n5-1,n6,topo)
		  + Sinv[k][6] * Int(d-2,n1,n2,n3,n4,n5,n6-1,topo) ) 
		 / ((d-sigma+n6-1) * b[6]) );


  }
 
  if (2 == d/2 - sigma + n ) {
    /*
     * Eq. 2.18
     */ 
    if ( n1 > 1 ) 
      RETURN6(( - Sinv[1][1] * Int(d-2,n1-1-1,n2,n3,n4,n5,n6,topo)
		- Sinv[1][2] * Int(d-2,n1-1,n2-1,n3,n4,n5,n6,topo)
		- Sinv[1][3] * Int(d-2,n1-1,n2,n3-1,n4,n5,n6,topo)
		- Sinv[1][4] * Int(d-2,n1-1,n2,n3,n4-1,n5,n6,topo)
		- Sinv[1][5] * Int(d-2,n1-1,n2,n3,n4,n5-1,n6,topo)
		- Sinv[1][6] * Int(d-2,n1-1,n2,n3,n4,n5,n6-1,topo)
		- b[1]*(d-sigma) * Int(d,n1-1,n2,n3,n4,n5,n6,topo) ) 
	      / (n1-1.0));
    
    if ( n2 > 1 ) 
      RETURN6(( - Sinv[2][1] * Int(d-2,n1-1,n2-1,n3,n4,n5,n6,topo)
		- Sinv[2][2] * Int(d-2,n1,n2-1-1,n3,n4,n5,n6,topo)
		- Sinv[2][3] * Int(d-2,n1,n2-1,n3-1,n4,n5,n6,topo)
		- Sinv[2][4] * Int(d-2,n1,n2-1,n3,n4-1,n5,n6,topo)
		- Sinv[2][5] * Int(d-2,n1,n2-1,n3,n4,n5-1,n6,topo)
		- Sinv[2][6] * Int(d-2,n1,n2-1,n3,n4,n5,n6-1,topo)
		- b[2]*(d-sigma) * Int(d,n1,n2-1,n3,n4,n5,n6,topo) ) 
	      / (n2-1.0));
    
    if ( n3 > 1 ) 
      RETURN6(( - Sinv[3][1] * Int(d-2,n1-1,n2,n3-1,n4,n5,n6,topo)
		- Sinv[3][2] * Int(d-2,n1,n2-1,n3-1,n4,n5,n6,topo)
		- Sinv[3][3] * Int(d-2,n1,n2,n3-1-1,n4,n5,n6,topo)
		- Sinv[3][4] * Int(d-2,n1,n2,n3-1,n4-1,n5,n6,topo)
		- Sinv[3][5] * Int(d-2,n1,n2,n3-1,n4,n5-1,n6,topo)
		- Sinv[3][6] * Int(d-2,n1,n2,n3-1,n4,n5,n6-1,topo)
		- b[3]*(d-sigma) * Int(d,n1,n2,n3-1,n4,n5,n6,topo) ) 
	      / (n3-1.0));
    
    if ( n4 > 1 ) 
      RETURN6(( - Sinv[4][1] * Int(d-2,n1-1,n2,n3,n4-1,n5,n6,topo)
		- Sinv[4][2] * Int(d-2,n1,n2-1,n3,n4-1,n5,n6,topo)
		- Sinv[4][3] * Int(d-2,n1,n2,n3-1,n4-1,n5,n6,topo)
		- Sinv[4][4] * Int(d-2,n1,n2,n3,n4-1-1,n5,n6,topo)
		- Sinv[4][5] * Int(d-2,n1,n2,n3,n4-1,n5-1,n6,topo)
		- Sinv[4][6] * Int(d-2,n1,n2,n3,n4-1,n5,n6-1,topo)
		- b[4]*(d-sigma) * Int(d,n1,n2,n3,n4-1,n5,n6,topo) ) 
	      / (n4-1.0));
    
    if ( n5 > 1 ) 
      RETURN6(( - Sinv[5][1] * Int(d-2,n1-1,n2,n3,n4,n5-1,n6,topo)
		- Sinv[5][2] * Int(d-2,n1,n2-1,n3,n4,n5-1,n6,topo)
		- Sinv[5][3] * Int(d-2,n1,n2,n3-1,n4,n5-1,n6,topo)
		- Sinv[5][4] * Int(d-2,n1,n2,n3,n4-1,n5-1,n6,topo)
		- Sinv[5][5] * Int(d-2,n1,n2,n3,n4,n5-1-1,n6,topo)
		- Sinv[5][5] * Int(d-2,n1,n2,n3,n4,n5-1,n6-1,topo)
		- b[5]*(d-sigma) * Int(d,n1,n2,n3,n4,n5-1,n6,topo) ) 
	      / (n5-1.0));

    if ( n6 > 1 ) 
      RETURN6(( - Sinv[5][1] * Int(d-2,n1-1,n2,n3,n4,n5,n6-1,topo)
		- Sinv[5][2] * Int(d-2,n1,n2-1,n3,n4,n5,n6-1,topo)
		- Sinv[5][3] * Int(d-2,n1,n2,n3-1,n4,n5,n6-1,topo)
		- Sinv[5][4] * Int(d-2,n1,n2,n3,n4-1,n5,n6-1,topo)
		- Sinv[5][5] * Int(d-2,n1,n2,n3,n4,n5-1,n6-1,topo)
		- Sinv[5][5] * Int(d-2,n1,n2,n3,n4,n5,n6-1-1,topo)
		- b[6]*(d-sigma) * Int(d,n1,n2,n3,n4,n5,n6-1,topo) ) 
	      / (n6-1.0));


  }

  cout << "# We should never end up here..." << endl;
  cout << "# Error in reduction of 6-point integrals." << endl;
  cout << "# Cannot reduce I(" << d << "," << n1 << "," << n2 << "," << n3 
       << "," << n4 << "," << n5 << "," << n6 << ")." << endl;
  exit(1);
}
#endif

/*********************************************************************
 * FIVE-POINT REDUCTION
 *********************************************************************/
IntType Int( int d,  int n1,  int n2,  int n3, 
	     int n4,  int n5, Topology* topo){


  if ( topo == 0 ) {
    cout << "# GGReduction: Null topology encountered" << endl;
    /***
     *** For testing purposes: return zero for undefined topologies.
     ***/ 
    //    return(0.);
    exit(1);
  }

#ifdef GGDEBUG
  cout << "# Pentagon(" << d << "," << n1 << "," << n2 << ","<< n3 << "," 
       << n4 << "," << n5 <<")" << endl;
#endif

#ifdef WITHCACHE

  /*
   * Check in the integral cache if we know this integral already:
   */

  map<unsigned int,IntType>::iterator integral = 
    (topo->cache).find(ihash(d,n1,n2,n3,n4,n5));
  
  if ( integral != (topo->cache).end() ) {
    return( (*integral).second );
  }

#endif

#ifdef NEWREC
  //--------------------------------------------------------------
  if (n1==0) RETURN5(Int(d,n2,n3,n4,n5,topo->subtopo(1)));
  if (n2==0) RETURN5(Int(d,n1,n3,n4,n5,topo->subtopo(2)));
  if (n3==0) RETURN5(Int(d,n1,n2,n4,n5,topo->subtopo(3)));
  if (n4==0) RETURN5(Int(d,n1,n2,n3,n5,topo->subtopo(4)));
  if (n5==0) RETURN5(Int(d,n1,n2,n3,n4,topo->subtopo(5)));
  //--------------------------------------------------------------
#endif

  const  int n = 5;
  const  int sigma = n1 + n2 + n3 + n4 + n5;

  double B, b[n+1], Sinv[n+1][n+1];

  B = topo->getB();


  if ( fabs(B) < BCUT ) {
#ifdef NEWREC
    double tmp = XInt(d,n1,n2,n3,n4,n5,topo,0);
    cout << "# new rec = " << tmp << endl;
    return(tmp);
#endif
    // We need a special reduction to avoid numerical instabilities:
#ifdef WITHEXCEPTIONS
    ostringstream message;
    message << "# 5-point reduction: \n";
    message << "# Special reduction not yet implemented. B-cut = ";
    message << BCUT;
    throw( GGException( message.str() ) );
#else
    cout << "# 5-point reduction:" << endl;
    cout << "# B = " << B << endl;
    cout << "# Special reduction not yet implemented..." <<endl;
#endif
    return 0.0;
  }
  
  
  if (n1==0) RETURN5(Int(d,n2,n3,n4,n5,topo->subtopo(1)));
  if (n2==0) RETURN5(Int(d,n1,n3,n4,n5,topo->subtopo(2)));
  if (n3==0) RETURN5(Int(d,n1,n2,n4,n5,topo->subtopo(3)));
  if (n4==0) RETURN5(Int(d,n1,n2,n3,n5,topo->subtopo(4)));
  if (n5==0) RETURN5(Int(d,n1,n2,n3,n4,topo->subtopo(5)));


  if ( d==6 && n1==1 && n2==1 && n3==1 && n4==1 && n5==1) {
    /***
     *** Evaluate pentagon integral in 6 dimensions, the final result should
     *** never depend on this value so we can set to whatever value we like:
     ***/
    RETURN5(0.0);
  }



  for (int i=1; i<= n; i++) {
    b[i] = 0.;
    for (int j=1; j<= n; j++) {
      Sinv[i][j] = topo->getSinv(i,j);
      b[i] += Sinv[i][j];
    }
  }



  /***
   ***  Magenta:
   ***/
  //  if (d==4 && n1>0 && n2>0 && n3>0 && n4>0 && n5>0) {
  if ( d==4 && n1==1 && n2==1 && n3==1 && n4==1 && n5==1 ) {
    RETURN5( ( 4.0 + 1.0 - sigma ) * B * Int(6,n1,n2,n3,n4,n5,topo)
	    + b[1]*Int(4,n1-1,n2,n3,n4,n5,topo)
	    + b[2]*Int(4,n1,n2-1,n3,n4,n5,topo)
	    + b[3]*Int(4,n1,n2,n3-1,n4,n5,topo)
	    + b[4]*Int(4,n1,n2,n3,n4-1,n5,topo)
	    + b[5]*Int(4,n1,n2,n3,n4,n5-1,topo));
  }

  
  /***
   ***  Red:
   ***/
  if ( d/2 + 3  == sigma ) {

    if ( n1 > 1 ) 
      RETURN5(( - Sinv[1][1] * Int(d-2,n1-1-1,n2,n3,n4,n5,topo)
		- Sinv[1][2] * Int(d-2,n1-1,n2-1,n3,n4,n5,topo)
		- Sinv[1][3] * Int(d-2,n1-1,n2,n3-1,n4,n5,topo)
		- Sinv[1][4] * Int(d-2,n1-1,n2,n3,n4-1,n5,topo)
		- Sinv[1][5] * Int(d-2,n1-1,n2,n3,n4,n5-1,topo)
		- b[1]*(d-sigma) * Int(d,n1-1,n2,n3,n4,n5,topo) ) / (n1-1.0));
    
    if ( n2 > 1 ) 
      RETURN5(( - Sinv[2][1] * Int(d-2,n1-1,n2-1,n3,n4,n5,topo)
		- Sinv[2][2] * Int(d-2,n1,n2-1-1,n3,n4,n5,topo)
		- Sinv[2][3] * Int(d-2,n1,n2-1,n3-1,n4,n5,topo)
		- Sinv[2][4] * Int(d-2,n1,n2-1,n3,n4-1,n5,topo)
		- Sinv[2][5] * Int(d-2,n1,n2-1,n3,n4,n5-1,topo)
		- b[2]*(d-sigma) * Int(d,n1,n2-1,n3,n4,n5,topo) ) / (n2-1.0));
    
    if ( n3 > 1 ) 
      RETURN5(( - Sinv[3][1] * Int(d-2,n1-1,n2,n3-1,n4,n5,topo)
		- Sinv[3][2] * Int(d-2,n1,n2-1,n3-1,n4,n5,topo)
		- Sinv[3][3] * Int(d-2,n1,n2,n3-1-1,n4,n5,topo)
		- Sinv[3][4] * Int(d-2,n1,n2,n3-1,n4-1,n5,topo)
		- Sinv[3][5] * Int(d-2,n1,n2,n3-1,n4,n5-1,topo)
		- b[3]*(d-sigma) * Int(d,n1,n2,n3-1,n4,n5,topo) ) / (n3-1.0));
    
    if ( n4 > 1 ) 
      RETURN5(( - Sinv[4][1] * Int(d-2,n1-1,n2,n3,n4-1,n5,topo)
		- Sinv[4][2] * Int(d-2,n1,n2-1,n3,n4-1,n5,topo)
		- Sinv[4][3] * Int(d-2,n1,n2,n3-1,n4-1,n5,topo)
		- Sinv[4][4] * Int(d-2,n1,n2,n3,n4-1-1,n5,topo)
		- Sinv[4][5] * Int(d-2,n1,n2,n3,n4-1,n5-1,topo)
		- b[4]*(d-sigma) * Int(d,n1,n2,n3,n4-1,n5,topo) ) / (n4-1.0));
    
    if ( n5 > 1 ) 
      RETURN5(( - Sinv[5][1] * Int(d-2,n1-1,n2,n3,n4,n5-1,topo)
		- Sinv[5][2] * Int(d-2,n1,n2-1,n3,n4,n5-1,topo)
		- Sinv[5][3] * Int(d-2,n1,n2,n3-1,n4,n5-1,topo)
		- Sinv[5][4] * Int(d-2,n1,n2,n3,n4-1,n5-1,topo)
		- Sinv[5][5] * Int(d-2,n1,n2,n3,n4,n5-1-1,topo)
		- b[5]*(d-sigma) * Int(d,n1,n2,n3,n4,n5-1,topo) ) / (n5-1.0));
  }

  /***
   ***  Blue:
   ***/
  if ( ( d/2 + 1  == sigma ) && (  d != sigma + 1 ) ) {
    RETURN5( ( Int(d-2,n1,n2,n3,n4,n5,topo)
	      - b[1]*Int(d-2,n1-1,n2,n3,n4,n5,topo)
	      - b[2]*Int(d-2,n1,n2-1,n3,n4,n5,topo)
	      - b[3]*Int(d-2,n1,n2,n3-1,n4,n5,topo)
	      - b[4]*Int(d-2,n1,n2,n3,n4-1,n5,topo)
	      - b[5]*Int(d-2,n1,n2,n3,n4,n5-1,topo) )/B/(d-1.0-sigma));
  }

  /***
   ***  Green:
   ***/
  if ( d/2 + 2  == sigma ) {
    if ( n1 > 1 ) {
      RETURN5( (-b[1]/B*Int(d-2,n1-1,n2,n3,n4,n5,topo) 
	       + (b[1]*b[1]/B - Sinv[1][1]) * Int(d-2,n1-2,n2,n3,n4,n5,topo) 
	       + (b[1]*b[2]/B - Sinv[1][2]) * Int(d-2,n1-1,n2-1,n3,n4,n5,topo) 
	       + (b[1]*b[3]/B - Sinv[1][3]) * Int(d-2,n1-1,n2,n3-1,n4,n5,topo)
	       + (b[1]*b[4]/B - Sinv[1][4]) * Int(d-2,n1-1,n2,n3,n4-1,n5,topo)
	       + (b[1]*b[5]/B - Sinv[1][5]) * Int(d-2,n1-1,n2,n3,n4,n5-1,topo)
	       )/(n1-1.0));
    };

    if ( n2 > 1 ) {
      RETURN5((-b[2]/B*Int(d-2,n1,n2-1,n3,n4,n5,topo) 
	      + (b[2]*b[1]/B - Sinv[2][1]) * Int(d-2,n1-1,n2-1,n3,n4,n5,topo) 
	      + (b[2]*b[2]/B - Sinv[2][2]) * Int(d-2,n1,n2-2,n3,n4,n5,topo) 
	      + (b[2]*b[3]/B - Sinv[2][3]) * Int(d-2,n1,n2-1,n3-1,n4,n5,topo)
	      + (b[2]*b[4]/B - Sinv[2][4]) * Int(d-2,n1,n2-1,n3,n4-1,n5,topo)
	      + (b[2]*b[5]/B - Sinv[2][5]) * Int(d-2,n1,n2-1,n3,n4,n5-1,topo)
	      )/(n2-1.0));
    };

    if ( n3 > 1 ) {
      RETURN5( (-b[3]/B*Int(d-2,n1,n2,n3-1,n4,n5,topo) 
	       + (b[3]*b[1]/B - Sinv[3][1]) * Int(d-2,n1-1,n2,n3-1,n4,n5,topo) 
	       + (b[3]*b[2]/B - Sinv[3][2]) * Int(d-2,n1,n2-1,n3-1,n4,n5,topo) 
	       + (b[3]*b[3]/B - Sinv[3][3]) * Int(d-2,n1,n2,n3-2,n4,n5,topo)
	       + (b[3]*b[4]/B - Sinv[3][4]) * Int(d-2,n1,n2,n3-1,n4-1,n5,topo)
	       + (b[3]*b[5]/B - Sinv[3][5]) * Int(d-2,n1,n2,n3-1,n4,n5-1,topo)
	       )/(n3-1.0));
    };

    if ( n4 > 1 ) {
      RETURN5((-b[4]/B*Int(d-2,n1,n2,n3,n4-1,n5,topo) 
	      + (b[4]*b[1]/B - Sinv[4][1]) * Int(d-2,n1-1,n2,n3,n4-1,n5,topo) 
	      + (b[4]*b[2]/B - Sinv[4][2]) * Int(d-2,n1,n2-1,n3,n4-1,n5,topo) 
	      + (b[4]*b[3]/B - Sinv[4][3]) * Int(d-2,n1,n2,n3-1,n4-1,n5,topo)
	      + (b[4]*b[4]/B - Sinv[4][4]) * Int(d-2,n1,n2,n3,n4-2,n5,topo)
	      + (b[4]*b[5]/B - Sinv[4][5]) * Int(d-2,n1,n2,n3,n4-1,n5-1,topo)
	      )/(n4-1.0));
    };

    if ( n5 > 1 ) {
      RETURN5((-b[5]/B*Int(d-2,n1,n2,n3,n4,n5-1,topo) 
	      + (b[5]*b[1]/B - Sinv[5][1]) * Int(d-2,n1-1,n2,n3,n4,n5-1,topo) 
	      + (b[5]*b[2]/B - Sinv[5][2]) * Int(d-2,n1,n2-1,n3,n4,n5-1,topo) 
	      + (b[5]*b[3]/B - Sinv[5][3]) * Int(d-2,n1,n2,n3-1,n4,n5-1,topo)
	      + (b[5]*b[4]/B - Sinv[5][4]) * Int(d-2,n1,n2,n3,n4-1,n5-1,topo)
	      + (b[5]*b[5]/B - Sinv[5][5]) * Int(d-2,n1,n2,n3,n4,n5-2,topo)
	      )/(n5-1.0));
    };
  }

  /*
   * Somehow I(10,1,1,1,1,1) was not reduced when applying the library
   * to the Hexagons, so I added the following.
   */
  if ( (  d != sigma + 1 ) ) {
    RETURN5( ( Int(d-2,n1,n2,n3,n4,n5,topo)
	      - b[1]*Int(d-2,n1-1,n2,n3,n4,n5,topo)
	      - b[2]*Int(d-2,n1,n2-1,n3,n4,n5,topo)
	      - b[3]*Int(d-2,n1,n2,n3-1,n4,n5,topo)
	      - b[4]*Int(d-2,n1,n2,n3,n4-1,n5,topo)
	      - b[5]*Int(d-2,n1,n2,n3,n4,n5-1,topo) )/B/(d-1.0-sigma));
  }

  cout << "# We should never end up here..." << endl;
  cout << "# Error in reduction of 5-point integrals." << endl;
  cout << "# Cannot reduce I(" << d << "," << n1 << "," << n2 << "," << n3 
       << "," << n4 << "," << n5 << ")." << endl;
  exit(1);
}

/*********************************************************************
 * FOUR-POINT REDUCTION
 *********************************************************************/
IntType Int( int d,  int n1,  int n2, 
	     int n3,  int n4, Topology* topo){


  if ( topo == 0 ) {
    cout << "# GGReduction: Null topology encountered" << endl;
    //    return(0.);
    exit(1);
  }


#ifdef WITHCACHE

  map<unsigned int,IntType>::iterator integral = 
    (topo->cache).find(ihash(d,n1,n2,n3,n4));
  
  if ( integral != (topo->cache).end() ) {
    return( (*integral).second );
  }

#endif

#ifdef GGDEBUG
  cout << "# Box(" << d << "," << n1 << "," << n2 << "," << n3 << "," 
       << n4 <<")" << endl;
#endif

  const  int n=4;
  const  int sigma = n1 + n2 + n3 + n4;

  double B,b[n+1],Sinv[n+1][n+1];


  B = topo->getB();
 
  //  cout << "# 4point: det(S) = " << topo->getDetS() << endl;

  if ( fabs(B) < BCUT ) {
    // Special reduction needed to avoid numerical problems
#ifdef WITHEXCEPTIONS
    ostringstream message;
    message << "# 4-point reduction: ";
    message << "# Special reduction not yet implemented. B-cut = ";
    message << BCUT;
    throw( GGException( message.str() ) );
#else
    cout << "# 4-point reduction:" << endl;
    cout << "# B = " << B << endl;
    cout << "# Special reduction not yet implemented..."<<endl;
#endif
  }



  if ( n1 == 0 ) RETURN4(Int(d,n2,n3,n4,topo->subtopo(1)));
  if ( n2 == 0 ) RETURN4(Int(d,n1,n3,n4,topo->subtopo(2)));
  if ( n3 == 0 ) RETURN4(Int(d,n1,n2,n4,topo->subtopo(3)));
  if ( n4 == 0 ) RETURN4(Int(d,n1,n2,n3,topo->subtopo(4)));


  for ( int i=1; i<= n; i++) {
    b[i] = 0.;
    for ( int j=1; j<= n; j++) {
      Sinv[i][j] = topo->getSinv(i,j);
      b[i] += Sinv[i][j];
    }
  }


  //  if ( d==4 && n1==1 && n2==1 && n3==1 && n4==1 ) {
  //  RETURN4((ScalarInt(topo->getMomenta(),topo->getMasses())));
  //  }

  if ( d==6 && n1==1 && n2==1 && n3==1 && n4==1 ) {
    
    RETURN4(
	    (ScalarInt(Dinput( 
			      topo->getS(1,2) + topo->getMass2(1) + topo->getMass2(2),
			      topo->getS(2,3) + topo->getMass2(2) + topo->getMass2(3), 
			      topo->getS(3,4) + topo->getMass2(3) + topo->getMass2(4),
			      topo->getS(1,4) + topo->getMass2(1) + topo->getMass2(4),
			      topo->getS(1,3) + topo->getMass2(1) + topo->getMass2(3),
			      topo->getS(2,4) + topo->getMass2(2) + topo->getMass2(4),
			      topo->getMass2(1), topo->getMass2(2),
			      topo->getMass2(3),topo->getMass2(4))
		       )
	     -b[1]*Int(4,1,1,1,topo->subtopo(1))
	     -b[2]*Int(4,1,1,1,topo->subtopo(2))
	     -b[3]*Int(4,1,1,1,topo->subtopo(3))
	     -b[4]*Int(4,1,1,1,topo->subtopo(4))
	    )/B/(d-1.0-sigma) );

    //return(-(b[1]+b[2]+b[3]+b[4])/B);
    //return(0.0); 
    //    return(1.0/B); 
  }
  

  /*
   *  Magenta:
   */
  if ( d == 4 && n1 > 0 && n2 > 0 && n3 > 0 && n4 > 0 ) {

#ifdef GGDEBUG    
    cout << "# 4-point Magenta\n";
#endif

    RETURN4( (d+1.0-sigma)*B*Int(6,n1,n2,n3,n4,topo)
	     + b[1]*Int(4,n1-1,n2,n3,n4,topo)
	     + b[2]*Int(4,n1,n2-1,n3,n4,topo)
	     + b[3]*Int(4,n1,n2,n3-1,n4,topo)
	     + b[4]*Int(4,n1,n2,n3,n4-1,topo)
	     );
  }
  
  /*
   *  Red:
   */
  if ( d/2 + 2  == sigma ) {

#ifdef GGDEBUG    
    cout << "# 4-point Red\n";
#endif

    if ( n1 > 1 ) 
      RETURN4((
	       -Sinv[1][1] * Int(d-2,n1-1-1,n2,n3,n4,topo)
	       -Sinv[1][2] * Int(d-2,n1-1,n2-1,n3,n4,topo)
	       -Sinv[1][3] * Int(d-2,n1-1,n2,n3-1,n4,topo)
	       -Sinv[1][4] * Int(d-2,n1-1,n2,n3,n4-1,topo)
	       - b[1]*(d-sigma)*Int(d,n1-1,n2,n3,n4,topo) ) / (n1-1.0));
    
    if ( n2 > 1 ) 
      RETURN4((
	       -Sinv[2][1] * Int(d-2,n1-1,n2-1,n3,n4,topo)
	       -Sinv[2][2] * Int(d-2,n1,n2-1-1,n3,n4,topo)
	       -Sinv[2][3] * Int(d-2,n1,n2-1,n3-1,n4,topo)
	       -Sinv[2][4] * Int(d-2,n1,n2-1,n3,n4-1,topo)
	       - b[2]*(d-sigma)*Int(d,n1,n2-1,n3,n4,topo) ) / (n2-1.0));
    
    if ( n3 > 1 ) 
      RETURN4((
	       -Sinv[3][1] * Int(d-2,n1-1,n2,n3-1,n4,topo)
	       -Sinv[3][2] * Int(d-2,n1,n2-1,n3-1,n4,topo)
	       -Sinv[3][3] * Int(d-2,n1,n2,n3-1-1,n4,topo)
	       -Sinv[3][4] * Int(d-2,n1,n2,n3-1,n4-1,topo)
	       - b[3]*(d-sigma)*Int(d,n1,n2,n3-1,n4,topo) ) / (n3-1.0));
    
    if ( n4 > 1 ) 
      RETURN4(( 
	       -Sinv[4][1] * Int(d-2,n1-1,n2,n3,n4-1,topo)
	       -Sinv[4][2] * Int(d-2,n1,n2-1,n3,n4-1,topo)
	       -Sinv[4][3] * Int(d-2,n1,n2,n3-1,n4-1,topo)
	       -Sinv[4][4] * Int(d-2,n1,n2,n3,n4-1-1,topo)
	       - b[4]*(d-sigma)*Int(d,n1,n2,n3,n4-1,topo)) / (n4-1.0));
  }

  /*
   *  Blue:
   */
  if ( ( d/2 == sigma ) && ( d != sigma + 1 ) ) {

#ifdef GGDEBUG    
    cout << "# 4-point Blue\n";
#endif

    double UVfinite = -2.0/(d-1.0-sigma)/(d-1.0-sigma) 
      * pow(-1.0,sigma-1)/gam[sigma-1];

    RETURN4(( Int(d-2,n1,n2,n3,n4,topo)
	      - b[1]*Int(d-2,n1-1,n2,n3,n4,topo)
	      - b[2]*Int(d-2,n1,n2-1,n3,n4,topo)
	      - b[3]*Int(d-2,n1,n2,n3-1,n4,topo)
	      - b[4]*Int(d-2,n1,n2,n3,n4-1,topo)
	     )/B/(d-1.0-sigma) + UVfinite*expflag);
  }

  /*
   *  Green:
   */
  if ( d/2+1  == sigma ) {

#ifdef GGDEBUG    
    cout << "# 4-point Green";
#endif
    if ( n1 > 1 ) {

      RETURN4( (-b[1]/B*Int(d-2,n1-1,n2,n3,n4,topo) 
		+ (b[1]*b[1]/B - Sinv[1][1]) * (Int(d-2,n1-2,n2,n3,n4,topo)
						-Int(d-2,n1-1,n2,n3,n4-1,topo)) 
		+ (b[1]*b[2]/B - Sinv[1][2]) * (Int(d-2,n1-1,n2-1,n3,n4,topo)
						-Int(d-2,n1-1,n2,n3,n4-1,topo)) 
		+ (b[1]*b[3]/B - Sinv[1][3]) * (Int(d-2,n1-1,n2,n3-1,n4,topo)
						-Int(d-2,n1-1,n2,n3,n4-1,topo))
	       )/(n1-1.0));

      /* The version above should be numerically more stable
      RETURN4( (-b[1]/B*Int(d-2,n1-1,n2,n3,n4,topo) 
		+ (b[1]*b[1]/B - Sinv[1][1]) * Int(d-2,n1-2,n2,n3,n4,topo) 
		+ (b[1]*b[2]/B - Sinv[1][2]) * Int(d-2,n1-1,n2-1,n3,n4,topo) 
		+ (b[1]*b[3]/B - Sinv[1][3]) * Int(d-2,n1-1,n2,n3-1,n4,topo)
		+ (b[1]*b[4]/B - Sinv[1][4]) * Int(d-2,n1-1,n2,n3,n4-1,topo)
		)/(n1-1.0));
      */
    };

    if ( n2 > 1 ) {

      RETURN4( (-b[2]/B*Int(d-2,n1,n2-1,n3,n4,topo) 
		+ (b[2]*b[1]/B - Sinv[2][1]) * (Int(d-2,n1-1,n2-1,n3,n4,topo)
						-Int(d-2,n1,n2-1,n3,n4-1,topo)) 
		+ (b[2]*b[2]/B - Sinv[2][2]) * (Int(d-2,n1,n2-2,n3,n4,topo)
						-Int(d-2,n1,n2-1,n3,n4-1,topo)) 
		+ (b[2]*b[3]/B - Sinv[2][3]) * (Int(d-2,n1,n2-1,n3-1,n4,topo)
						-Int(d-2,n1,n2-1,n3,n4-1,topo))
		)/(n2-1.0));

      /* The version above should be numerically more stable
      RETURN4( (-b[2]/B*Int(d-2,n1,n2-1,n3,n4,topo) 
		+ (b[2]*b[1]/B - Sinv[2][1]) * Int(d-2,n1-1,n2-1,n3,n4,topo) 
		+ (b[2]*b[2]/B - Sinv[2][2]) * Int(d-2,n1,n2-2,n3,n4,topo) 
		+ (b[2]*b[3]/B - Sinv[2][3]) * Int(d-2,n1,n2-1,n3-1,n4,topo)
		+ (b[2]*b[4]/B - Sinv[2][4]) * Int(d-2,n1,n2-1,n3,n4-1,topo)
		)/(n2-1.0));
      */
    };

    if ( n3 > 1 ) {
      RETURN4( (- b[3]/B*Int(d-2,n1,n2,n3-1,n4,topo) 
		+ (b[3]*b[1]/B - Sinv[3][1]) * (Int(d-2,n1-1,n2,n3-1,n4,topo)
						-Int(d-2,n1,n2,n3-1,n4-1,topo)) 
		+ (b[3]*b[2]/B - Sinv[3][2]) * (Int(d-2,n1,n2-1,n3-1,n4,topo)
						-Int(d-2,n1,n2,n3-1,n4-1,topo)) 
		+ (b[3]*b[3]/B - Sinv[3][3]) * (Int(d-2,n1,n2,n3-2,n4,topo)
						-Int(d-2,n1,n2,n3-1,n4-1,topo))
		)/(n3-1.0));

      /* The version above should be numerically more stable
      RETURN4( (-b[3]/B*Int(d-2,n1,n2,n3-1,n4,topo) 
		+ (b[3]*b[1]/B - Sinv[3][1]) * Int(d-2,n1-1,n2,n3-1,n4,topo) 
		+ (b[3]*b[2]/B - Sinv[3][2]) * Int(d-2,n1,n2-1,n3-1,n4,topo) 
		+ (b[3]*b[3]/B - Sinv[3][3]) * Int(d-2,n1,n2,n3-2,n4,topo)
		+ (b[3]*b[4]/B - Sinv[3][4]) * Int(d-2,n1,n2,n3-1,n4-1,topo)
	       )/(n3-1.0));
      */
    };

    if ( n4 > 1 ) {

      RETURN4( (-b[4]/B*Int(d-2,n1,n2,n3,n4-1,topo) 
		+ (b[4]*b[1]/B - Sinv[4][1]) * (Int(d-2,n1-1,n2,n3,n4-1,topo)
						-Int(d-2,n1,n2,n3,n4-2,topo)) 
		+ (b[4]*b[2]/B - Sinv[4][2]) * (Int(d-2,n1,n2-1,n3,n4-1,topo)
						-Int(d-2,n1,n2,n3,n4-2,topo))  
		+ (b[4]*b[3]/B - Sinv[4][3]) * (Int(d-2,n1,n2,n3-1,n4-1,topo)
						-Int(d-2,n1,n2,n3,n4-2,topo)) 
		)/(n4-1.0));

      /* The version above should be numerically more stable

      RETURN4( (-b[4]/B*Int(d-2,n1,n2,n3,n4-1,topo) 
		+ (b[4]*b[1]/B - Sinv[4][1]) * Int(d-2,n1-1,n2,n3,n4-1,topo) 
		+ (b[4]*b[2]/B - Sinv[4][2]) * Int(d-2,n1,n2-1,n3,n4-1,topo) 
		+ (b[4]*b[3]/B - Sinv[4][3]) * Int(d-2,n1,n2,n3-1,n4-1,topo)
		+ (b[4]*b[4]/B - Sinv[4][4]) * Int(d-2,n1,n2,n3,n4-2,topo)
		)/(n4-1.0));
      */
    };
  }

#ifdef NEWREC
  // epsilon/epsilon terms need to be fixed
  return( ( + Int(d-2,n1,n2,n3,n4,topo) 
	    - b[1] * Int(d-2,n1-1,n2,n3,n4,topo)
	    - b[2] * Int(d-2,n1,n2-1,n3,n4,topo)
	    - b[3] * Int(d-2,n1,n2,n3-1,n4,topo)
	    - b[4] * Int(d-2,n1,n2,n3,n4-1,topo) 
	    ) / (d-1.-sigma) / B );
#endif

  cout << "# We should never end up here..." << endl;
  cout << "# Error in reduction of 4-point integrals." << endl;
  cout << "# Cannot reduce I(" << d << "," << n1 << "," << n2 << "," << n3 
       << "," << n4 << ")." << endl;
  exit(1);
}

/******************************************************************
 *
 * Reduction of the 3-point topologies
 *
 *
 ******************************************************************/
IntType Int( int d,  int n1,  int n2, 
	     int n3, Topology* topo){

#ifdef GGDEBUG
  cout << "# Triangle(" << d << "," << n1 << "," << n2 << ","<<n3<<")" << endl;
#endif

#ifdef WITHCACHE
  map<unsigned int,IntType>::iterator integral = 
    (topo->cache).find(ihash(d,n1,n2,n3));
  
  if ( integral != (topo->cache).end() ) {
    return( (*integral).second );
  }
#endif

  const  int n=3;
  const  int sigma=n1+n2+n3;
  double B,b[n+1],Sinv[n+1][n+1];

  if ( n1 == 0 ) {
    Topology* t2 = topo->subtopo(1);
    double m1q= t2->getMass2(1);
    double m2q= t2->getMass2(2);
    double pq = t2->getS(1,2)+m1q+m2q;
    RETURN3(Int(d,n2,n3,pq,m1q,m2q));
  }
  if ( n2 == 0 ) {
    Topology* t2 = topo->subtopo(2);
    double m1q= t2->getMass2(1);
    double m2q= t2->getMass2(2);
    double pq = t2->getS(1,2)+m1q+m2q;
    RETURN3(Int(d,n1,n3,pq,m1q,m2q));}
  
  if ( n3 == 0 ) {
    Topology* t2 = topo->subtopo(3);
    double m1q= t2->getMass2(1);
    double m2q= t2->getMass2(2);
    double pq = t2->getS(1,2)+m1q+m2q;
    RETURN3(Int(d,n1,n2,pq,m1q,m2q));
  }



  if ( d==4 && n1==1 && n2==1 && n3==1 ) {
    RETURN3(ScalarInt( Cinput( topo->getS(1,2) + topo->getMass2(1) + topo->getMass2(2),
			       topo->getS(2,3) + topo->getMass2(2) + topo->getMass2(3), 
			       topo->getS(1,3) + topo->getMass2(1) + topo->getMass2(3),
			       topo->getMass2(1), topo->getMass2(2),
			       topo->getMass2(3) ) 
		       ) 
	    ); 
  }
  
  if ( ( topo->getDetS() == 0.0 ) ) {

    // cout << "# Triangle(" << d << "," << n1 << "," << n2 << ","<<n3<<")" << endl;
    // cout << "# new rec " << endl;
    // cout << "# finite terms " << endl;

    double z1,z2,z3;
    z1 = topo->getZ(1);
    z2 = topo->getZ(2);
    z3 = topo->getZ(3);
    double z = z1 + z2 + z3;
    
    if ( d/2 == sigma ) {
      // cout << "# new rec a" << endl;
      double UVfinite = -2.0/(d-1.0-sigma)/(d-1.0-sigma)
	* pow(-1.0,sigma-1)/gam[sigma-1];
      
      RETURN3( - ( + z1*Int(d-2,n1-1,n2,n3,topo)
		   + z2*Int(d-2,n1,n2-1,n3,topo) 
		   + z3*Int(d-2,n1,n2,n3-1,topo) )/(d-1.0-sigma)/z
	       + UVfinite * expflag);
    }
    
    
    if ( ( d/2 + 1 ) == sigma ) {

      // cout << "# new rec b" << endl;
      for ( int i=1; i<= n; i++) {
	b[i] = 0.;
	for ( int j=1; j<= n; j++) {
	  Sinv[i][j] = topo->getSinv(i,j);
	  b[i] += Sinv[i][j];
	}
      }

#define SMALL 1.e-8     
      if ( fabs( topo->getRange(1) ) < SMALL ){

	if ( n1 > 1 ) {

	  RETURN3(- (+ (Sinv[1][1] - z1/z*b[1])*( Int(d-2,n1-2,n2,n3,topo)-Int(d-2,n1-1,n2,n3-1,topo))
		     + (Sinv[2][1] - z2/z*b[1])*(Int(d-2,n1-1,n2-1,n3,topo)-Int(d-2,n1-1,n2,n3-1,topo)) 
		     )/(n1-1.0)
		  );
      /* The version above should be numerically more stable
	  RETURN3(- (+ (Sinv[1][1] - z1/z*b[1])*Int(d-2,n1-2,n2,n3,topo) 
		     + (Sinv[2][1] - z2/z*b[1])*Int(d-2,n1-1,n2-1,n3,topo) 
		     + (Sinv[3][1] - z3/z*b[1])*Int(d-2,n1-1,n2,n3-1,topo) 
		     )/(n1-1.0)
		  );
      */
	}
      }

      if (fabs(topo->getRange(2)) < SMALL) {

	if ( n2 > 1 ) {
	  
	  RETURN3(- ( + (Sinv[1][2] - z1/z*b[2])*(Int(d-2,n1-1,n2-1,n3,topo)
						  -Int(d-2,n1,n2-1,n3-1,topo)) 
		      + (Sinv[2][2] - z2/z*b[2])*(Int(d-2,n1,n2-2,n3,topo)
						  -Int(d-2,n1,n2-1,n3-1,topo)) 
		      )/(n2-1.0)
		      );
	  /* The version above should be numerically more stable
	    RETURN3(- ( + (Sinv[1][2] - z1/z*b[2])*Int(d-2,n1-1,n2-1,n3,topo) 
	    + (Sinv[2][2] - z2/z*b[2])*Int(d-2,n1,n2-2,n3,topo) 
	    + (Sinv[3][2] - z3/z*b[2])*Int(d-2,n1,n2-1,n3-1,topo) 
	    )/(n2-1.0)
	    );
	  */
	}
      }

      if (fabs(topo->getRange(3)) < SMALL) {
	
	if ( n3 > 1 ) {
	  RETURN3(- ( + (Sinv[1][3] - z1/z*b[3])*(Int(d-2,n1-1,n2,n3-1,topo)-Int(d-2,n1,n2,n3-2,topo))
		      + (Sinv[2][3] - z2/z*b[3])*(Int(d-2,n1,n2-1,n3-1,topo)- Int(d-2,n1,n2,n3-2,topo))
		      )/(n3-1.0)
		 );

	  RETURN3(- ( + (Sinv[1][3] - z1/z*b[3])*Int(d-2,n1-1,n2,n3-1,topo) 
		      + (Sinv[2][3] - z2/z*b[3])*Int(d-2,n1,n2-1,n3-1,topo) 
		      + (Sinv[3][3] - z3/z*b[3])*Int(d-2,n1,n2,n3-2,topo) 
		      )/(n3-1.0)
		  );

	}
      }
      /*
      if (( fabs(topo->getRange(1))+fabs(topo->getRange(2))) < 1.e-11){

	RETURN3( ( - Int(d-2,n1,n2,n3-1,topo)
		  - double(n1) * Int(d,n1+1,n2,n3-1,topo)
		  - double(n2) * Int(d,n1,n2+1,n3-1,topo))/(n3-1.0) );
      }

      if ((fabs(topo->getRange(1))+fabs(topo->getRange(3))) < 1.e-11){

	RETURN3( ( - Int(d-2,n1,n2-1,n3,topo)
		  - double(n1) * Int(d,n1+1,n2-1,n3,topo)
		  - double(n3) * Int(d,n1,n2-1,n3+1,topo))/(n2-1.0) );

      }

      if ((fabs(topo->getRange(2))+fabs(topo->getRange(3))) < 1.e-11){

	RETURN3( ( - Int(d-2,n1-1,n2,n3,topo)
		  - double(n2) * Int(d,n1-1,n2+1,n3,topo)
		  - double(n3) * Int(d,n1-1,n2,n3+1,topo))/(n1-1.0) );

      }

      */
      //cout << "# Triangle(" << d << "," << n1 << "," << n2 << ","<<n3<<")"  << endl;


      if ( fabs(topo->getRange(3)) < SMALL){
	double a[4],y[4]={0.,0.,0.,0.},ysum;
	a[1] = z2; 
	a[2] = -z1; 
	a[3] = 0.0;  

	for( int i=1; i<4;i++)
	  for( int j=1; j<4;j++)
	    y[i] += Sinv[i][j]*a[j];
	      
	ysum = y[1]+y[2]+y[3];

	y[1] = y[1] - z1/z*ysum;       
	y[2] = y[2] - z2/z*ysum;       
	y[3] = y[3] - z3/z*ysum;       

	//	 	cout << "#1-2\n";
	// 	cout << a[1] << "# " << a[2] << " " << a[3] << endl;
	// 	cout << topo->getS(1,1)*y[1]+topo->getS(1,2)*y[2]+topo->getS(1,3)*y[3] << " ";
	// 	cout << topo->getS(2,1)*y[1]+topo->getS(2,2)*y[2]+topo->getS(2,3)*y[3] << " " ;
	// 	cout << topo->getS(3,1)*y[1]+topo->getS(3,2)*y[2]+topo->getS(3,3)*y[3] << endl;
	// 	cout << y[1]+y[2]+y[3] << endl;
	
	if ( n1 > 1 ) {
	  RETURN3(
		 (-y[1]*Int(d-2,n1-2,n2,n3,topo)
		  -y[2]*Int(d-2,n1-1,n2-1,n3,topo)
		  -y[3]*Int(d-2,n1-1,n2,n3-1,topo)
		  + a[2] * Int(d-2,n1-1,n2,n3,topo)
		  + a[2]*n3*Int(d,n1-1,n2,n3+1,topo))/(a[1]-a[2])/(n1-1.0)
		 );
	}
	
	if ( n2 > 1 ) {
	  RETURN3(
		 (-y[1]*Int(d-2,n1-1,n2-1,n3,topo)
		  -y[2]*Int(d-2,n1,n2-2,n3,topo)
		  -y[3]*Int(d-2,n1,n2-1,n3-1,topo)
		  + a[1] * Int(d-2,n1,n2-1,n3,topo)
		  + a[1] * n3 * Int(d,n1,n2-1,n3+1,topo))/(a[2]-a[1])/(n2-1.0)
		 );
	}
      }

      if ( fabs(topo->getRange(2)) < SMALL){
	double a[4],y[4]={0.,0.,0.,0.},ysum;
	a[1] = z3; 
	a[2] = 0.0; 
	a[3] = -z1;  

	for( int i=1; i<4;i++)
	  for( int j=1; j<4;j++)
	    y[i] += Sinv[i][j]*a[j];
	      
	ysum = y[1]+y[2]+y[3];
	
	y[1] = y[1] - z1/z*ysum;       
	y[2] = y[2] - z2/z*ysum;       
	y[3] = y[3] - z3/z*ysum;       
	//	 	cout << "#1-3\n";
	// 	cout << a[1] << " " << a[2] << " " << a[3] << endl;
	// 	cout << topo->getS(1,1)*y[1]+topo->getS(1,2)*y[2]+topo->getS(1,3)*y[3] << " ";
	// 	cout << topo->getS(2,1)*y[1]+topo->getS(2,2)*y[2]+topo->getS(2,3)*y[3] << " ";
	// 	cout << topo->getS(3,1)*y[1]+topo->getS(3,2)*y[2]+topo->getS(3,3)*y[3] << endl;
	// 	cout << y[1]+y[2]+y[3] << endl;

	if ( n1 > 1) {
	  RETURN3(
		  (-y[1]*Int(d-2,n1-2,n2,n3,topo)
		   -y[2]*Int(d-2,n1-1,n2-1,n3,topo)
		   -y[3]*Int(d-2,n1-1,n2,n3-1,topo)
		   + a[3]*Int(d-2,n1-1,n2,n3,topo)
		   + a[3]*n2*Int(d,n1-1,n2+1,n3,topo))/(a[1]-a[3])/(n1-1.0)
		  );
	}

	if ( n3 > 1 ) {
	  RETURN3(
		  (-y[1]*Int(d-2,n1-1,n2,n3-1,topo)
		   -y[2]*Int(d-2,n1,n2-1,n3-1,topo)
		   -y[3]*Int(d-2,n1,n2,n3-2,topo)
		   + a[1]*Int(d-2,n1,n2,n3-1,topo)
		   + a[1]*n2*Int(d,n1,n2+1,n3-1,topo))/(a[3]-a[1])/(n3-1.0)
		  );
	}
      }


      if ( fabs(topo->getRange(1)) < SMALL){

	double a[4],y[4]={0.,0.,0.,0.},ysum;
	a[1] = 0.; 
	a[2] = z3; 
	a[3] = -z2;  

	for( int i=1; i<4;i++)
	  for( int j=1; j<4;j++)
	    y[i] += Sinv[i][j]*a[j];
	
	ysum = y[1]+y[2]+y[3];

	y[1] = y[1] - z1/z*ysum;       
	y[2] = y[2] - z2/z*ysum;       
	y[3] = y[3] - z3/z*ysum;       
	//		cout << "2-3\n";
	// 	cout << a[1] << " " << a[2] << " " << a[3] << endl;
	// 	cout << topo->getS(1,1)*y[1]+topo->getS(1,2)*y[2]+topo->getS(1,3)*y[3] << " ";
	// 	cout << topo->getS(2,1)*y[1]+topo->getS(2,2)*y[2]+topo->getS(2,3)*y[3] << " ";
	// 	cout << topo->getS(3,1)*y[1]+topo->getS(3,2)*y[2]+topo->getS(3,3)*y[3] << endl;
	// 	cout << y[1]+y[2]+y[3] << endl;
	if ( n2 > 1 ) {
	  RETURN3(
		  (-y[1]*Int(d-2,n1-1,n2-1,n3,topo)
		   -y[2]*Int(d-2,n1,n2-2,n3,topo)
		   -y[3]*Int(d-2,n1,n2-1,n3-1,topo)
		   + a[3]*Int(d-2,n1,n2-1,n3,topo)
		   + a[3]*n1*Int(d,n1+1,n2-1,n3,topo))/(a[2]-a[3])/(n2-1.0)
		  );
	}
	
	if (n3>1) {
	  RETURN3(
		  (-y[1]*Int(d-2,n1-1,n2,n3-1,topo)
		   -y[2]*Int(d-2,n1,n2-1,n3-1,topo)
		   -y[3]*Int(d-2,n1,n2,n3-2,topo)
		   + a[2]*Int(d-2,n1,n2,n3-1,topo)
		   + a[2]*n1*Int(d,n1+1,n2,n3-1,topo))/(a[3]-a[2])/(n3-1.0)
		  );
	}
      }



    }
  }
  
  B = topo->getB();

  //cout << "# Triangle: "<< B <<endl;
  if ( fabs(B) < BCUT ) {
#ifdef NEWREC
    return(XInt(d,n1,n2,n3,topo,0));
#endif
    // Special reduction needed to avoid numerical problems
#ifdef WITHEXCEPTIONS
    ostringstream message;
    message << "# 3-point reduction: ";
    message << "# Special reduction not yet implemented. B-cut = ";
    message << BCUT;
    throw( GGException( message.str() ) );
#else
    cout << "# 3-point reduction:" << endl;
    cout << "# B = " << B << endl;
    topo->printS();
    cout << "# 3point: det(S) = " << topo->getDetS() << endl;
    //    cout << "# q_0 = " << topo->q[0] << endl;
    // cout << "# q_1 = " << topo->q[1] << endl;
    //cout << "# q_2 = " << topo->q[2] << endl;
    cout << "# z1 = " << topo->getZ(1) << endl;
    cout << "# z2 = " << topo->getZ(2) << endl;
    cout << "# z3 = " << topo->getZ(3) << endl;
    cout << "# R1 = " << topo->getRange(1) << endl;
    cout << "# R2 = " << topo->getRange(2) << endl;
    cout << "# R3 = " << topo->getRange(3) << endl;
    cout << "# Special reduction not yet implemented..."<<endl;
    cout << "# when trying to reduce I(" << d << ","
	 << n1 << ","
	 << n2 << ","
	 << n3 << ")" << endl;
    exit(1);
#endif
  }


  for ( int i=1; i<= n; i++) {
    b[i] = 0.;
    for ( int j=1; j<= n; j++) {
      Sinv[i][j] = topo->getSinv(i,j);
      b[i] += Sinv[i][j];
    }
  }




  /*
   *  Blue: (we are on the UV-Line)
   */
  if ( ( d/2 == sigma ) && ( d != sigma + 1 ) ) {
#ifdef GGDEBUG
    cout << "# 3-point blue" << endl;
#endif

    double UVfinite = -2.0/(d-1.0-sigma)/(d-1.0-sigma) 
      * pow(-1.0,sigma-1)/gam[sigma-1];

    RETURN3(( Int(d-2,n1,n2,n3,topo)
	     - b[1]*Int(d-2,n1-1,n2,n3,topo)
	     - b[2]*Int(d-2,n1,n2-1,n3,topo)
	     - b[3]*Int(d-2,n1,n2,n3-1,topo))/B/(d-1.0-sigma)
	    + UVfinite * expflag );
  }

  /*
   *  Green: (
   */
  if ( sigma == d/2 + 1 ) {
#ifdef GGDEBUG
    cout << "# 3-point green" << endl;
#endif

    if ( n1 > 1 ) {

      RETURN3( (- b[1]/B*Int(d-2,n1-1,n2,n3,topo) 
		+ (b[1]*b[1]/B - Sinv[1][1]) * (Int(d-2,n1-2,n2,n3,topo)
						- Int(d-2,n1-1,n2,n3-1,topo))
		+ (b[1]*b[2]/B - Sinv[1][2]) * (Int(d-2,n1-1,n2-1,n3,topo)
						- Int(d-2,n1-1,n2,n3-1,topo))
		)/(n1-1.0));
      /* The above version should be numerically more stable...
      RETURN3( (-b[1]/B*Int(d-2,n1-1,n2,n3,topo) 
	       + (b[1]*b[1]/B - Sinv[1][1]) * Int(d-2,n1-2,n2,n3,topo) 
	       + (b[1]*b[2]/B - Sinv[1][2]) * Int(d-2,n1-1,n2-1,n3,topo) 
	       + (b[1]*b[3]/B - Sinv[1][3]) * Int(d-2,n1-1,n2,n3-1,topo)
	       )/(n1-1.0));
      */
    };

    if ( n2 > 1 ) {
      RETURN3((- b[2]/B*Int(d-2,n1,n2-1,n3,topo) 
	       + (b[2]*b[1]/B - Sinv[2][1]) * (Int(d-2,n1-1,n2-1,n3,topo)
					       -Int(d-2,n1,n2-1,n3-1,topo)) 
	       + (b[2]*b[2]/B - Sinv[2][2]) * (Int(d-2,n1,n2-2,n3,topo)
					       -Int(d-2,n1,n2-1,n3-1,topo)) 
	       )/(n2-1.0));
      /* The above version should be numerically more stable...
      RETURN3((-b[2]/B*Int(d-2,n1,n2-1,n3,topo) 
	      + (b[2]*b[1]/B - Sinv[2][1]) * Int(d-2,n1-1,n2-1,n3,topo) 
	      + (b[2]*b[2]/B - Sinv[2][2]) * Int(d-2,n1,n2-2,n3,topo) 
	      + (b[2]*b[3]/B - Sinv[2][3]) * Int(d-2,n1,n2-1,n3-1,topo)
	      )/(n2-1.0));
      */
    };

    if (n3>1) {  

      RETURN3( ( - b[3]/B*Int(d-2,n1,n2,n3-1,topo) 
		 + (b[3]*b[1]/B - Sinv[3][1]) * (Int(d-2,n1-1,n2,n3-1,topo)
						 -Int(d-2,n1,n2,n3-2,topo)) 
		 + (b[3]*b[2]/B - Sinv[3][2]) * (Int(d-2,n1,n2-1,n3-1,topo) 
						 -Int(d-2,n1,n2,n3-2,topo))
		 )/(n3-1.0));
      /* The above version should be numerically more stable...
      RETURN3( ( -b[3]/B*Int(d-2,n1,n2,n3-1,topo) 
	       + (b[3]*b[1]/B - Sinv[3][1]) * Int(d-2,n1-1,n2,n3-1,topo) 
	       + (b[3]*b[2]/B - Sinv[3][2]) * Int(d-2,n1,n2-1,n3-1,topo) 
	       + (b[3]*b[3]/B - Sinv[3][3]) * Int(d-2,n1,n2,n3-2,topo)
	       )/(n3-1.0));
      */
    };

  }

#ifdef NEWREC
  // epsilon/epsilon terms need to be fixed
  return( ( + Int(d-2,n1,n2,n3,topo) 
	    - b[1] * Int(d-2,n1-1,n2,n3,topo)
	    - b[2] * Int(d-2,n1,n2-1,n3,topo)
	    - b[3] * Int(d-2,n1,n2,n3-1,topo)
	    ) / (d-1.-sigma) / B );
#endif

  cout << "# We should never end up here..." << endl;
  cout << "# Error in reduction of 3-point integrals," << endl;
  cout << "# Cannot reduce I(" << d << "," << n1 << "," << n2 << "," << n3 
       << ")." << endl;
  exit(1);
}

/******************************************************************
 *
 * Reduction of the 2-point topologies
 *
 *
 ******************************************************************/
IntType Int( int d,  int n1,  int n2, 
	    double pq, double  m1q, double m2q){

#ifdef GGDEBUG
  cout << "# Int(" << d << "," << n1 << "," << n2 << ")" << endl;
#endif

  if ( n1 == 0 ) 
    return Int(d,n2,m2q);

  if ( n2 == 0 ) 
    return Int(d,n1,m1q);

  if ( ( n1 > 1 ) && ( d > 4 ) )
    return( -( + Int(d-2,n1-1,n2,pq,m1q,m2q) 
	       + static_cast<double>(n2) * Int(d,n1-1,n2+1,pq,m1q,m2q) 
	      ) / ( n1 - 1.0 ) 
	    );

  if ( ( n1 == 1 ) && ( n2 > 1 ) && ( d > 4 ) ){
    if ( fabs(pq) > 1.e-7) {
      return(
	     ( Int(d-2,n1,n2-2,pq,m1q,m2q) - Int(d-2,n1-1,n2-1,pq,m1q,m2q) 
	      + (m2q-m1q-pq) * Int(d-2,n1,n2-1,pq,m1q,m2q) ) /(n2-1.0)/2.0/pq);
    } else {
      if ( m1q == m2q ) {
	if ( m1q!=0 ) {
	  return( Int(d,n2+1,m2q) );
	} else {
	  if ( d/2 - n1 - n2 == 0){ 
	    return(pow(-1.0,n1+n2)/gam[n1+n2]*ScalarInt(0.,0.,0.));
	  } else {
	    cout << "# Do not know what to do with this integral..." << endl;
	    exit(1);
	  }
	}
      } else {
	// Additional reduction added on 09/2009:
	if ( (n1==1) && (n2==2) && (d==6) ) {
	  return( ( Int(d,1,1,0,m1q,m2q) + m1q*Int(4,1,1,0,m1q,m2q) 
		    - expflag * (m2q+m1q)/2. )
		  /(m2q-m1q) );
	}
	if ( (n1==1) && (n2==3) && (d==8) ) {
	  return((2.*Int(8,1,1,0.,m1q,m2q) + 2.*m1q* Int(6,1,1,0.,m1q,m2q) 
		  + m1q*m1q*Int(4,1,1,0,m1q,m2q) 
		  + expflag*(-2*m1q*(m2q-m1q)/2. - (m1q*m1q+m1q*m2q+m2q*m2q) 
			     ) 
		  )/(2*pow(m2q-m1q,2) ) );
	}
	cout << "# Two-point reduction: p^2 = 0, m1 <> m2 not yet implemented" 
	     << endl;
	cout << "# Int(d="<< d <<",n1="<<n1 <<",n2="<<n2 
	     <<",pq="<<pq <<",m1q="<<m1q <<",m2q="<<m2q<<")" << endl; 
	exit(1);
      } 
    };
  };
  
  if (d==6 && (n1==1) && (n2==1)) { 
    if (fabs(pq) > 1.e-7) {
      double UVfinite = 1.0/(d-3)/(d-3)*(pq-3.0*m1q-3.0*m2q);

      return(
	     -1.0/(d-3.0)*
	     ( Int(d-2,1,m2q) + 2.0*m1q*Int(d-2,1,1,pq,m1q,m2q)
	       + (m1q-m2q+pq) * 0.5/pq 
	       * ( Int(d-2,1,m1q) - Int(d-2,1,m2q)
		   + (m2q-m1q-pq) * Int(d-2,1,1,pq,m1q,m2q) 
		   ) 
	       ) 
	     + UVfinite*expflag
	     );
    } else {
      if ( m1q == m2q ) {
	return(-Int(4,1,m1q));
      } else {
	// Additional reduction added on 09/2009:
	return((Int(d,1,m1q)-Int(d,1,m2q))/(m1q-m2q));
	cout << "# Two-point reduction: p² ~ 0, m1 <> m2 not yet implemented" 
	     << endl;
	exit(1);
      } 
    }
  }

  if ( ( d == 4 ) && ( n1 == 1 ) && ( n2 == 1 ) ) 
    return(ScalarInt(pq,m1q,m2q)); 

  // Additional reduction added on 09/2009:
  if ( ( d == 8 ) && ( n1 == 1 ) && ( n2 == 1 ) && (pq==0.) ) 
    return((Int(d,1,m1q)-Int(d,1,m2q))/(m1q-m2q));
  cout << "# We should never end up here..." << endl;
  cout << "# Error in reduction of two-point integrals." << endl;
  cout << "# when evaluating I(" << d << "," << n1 << "," << n2 << ")" << endl; 
  cout << "# m1q = " << m1q  << " m2q = " << m2q << " pq = " << pq << endl; 
  exit(1);
}


IntType Int( int d,  int n1, double mq){

  /*
   * Use explicit formulae to evaluate the tadpole integrals:
   * 
   * Note: Only the real part is evaluated!
   */

  if ( d/2  >= n1 ) {
    if ( mq != 0 ) {
      int h = d/2-n1;
      return( pow(-1.,n1) * pow(-mq,h) / gam[n1] / gam[h+1] 
	      * ( parms.getDeltaUV1() 
		  + log( parms.getInternalScaleSquared()/fabs(mq) )
		  + S1[h] ) );
    }
  } else {
    return( pow(-1.,n1) * gam[n1-d/2] / gam[n1] * pow(mq,d/2-n1) );
  }
						    

  /*
   * Reduce both d/2 and n1 by one unit
   */
  if ( ( d > 4 ) && ( n1 > 1 ) ) {
    return( - Int(d-2,n1-1,mq) / ( n1 - 1.0 ) );
  }


  /*
   * Reduce n1 by one unit
   */
  if ( ( d >= 4 ) && ( n1 > 1 ) ) {
    double UVfinite = 0.;
    if ( n1 < d/2 + 2 )
      UVfinite = pow(-1.,d/2+1) / gam[d/2+2-n1]/gam[n1-1] / ( n1 - 1.0 ) 
	* pow(mq, d/2 - n1 ) * expflag;

    return( - ( n1 - 1.0 - d/2. ) / mq / ( n1 - 1.0 ) * Int(d,n1-1,mq) 
	   + UVfinite );
  }

  /*
   * Reduce both d/2 by one unit
   */
  if ( ( d > 4 ) && ( n1 >= 1 ) && ( n1 != d/2) ) {
    double UVfinite = 0.;
    if ( n1   < d/2 )
      UVfinite = 2.0 / ( 2.0*n1 - d ) * pow(-1.,d/2+1) 
	* pow(mq, d/2 - n1)/gam[d/2-n1+1]/gam[n1] * expflag;

    return( 1.0 / ( n1 - d/2.0 ) * mq * Int(d-2,n1,mq) + UVfinite );
  }
 

  if ( ( d == 4 ) && ( n1 == 1 ) ) 
    return(ScalarInt(mq));

  cout << "# Reduction of one-point integrals failed " << endl;
  cout << "# I(" << d << "," <<n1<<","<<mq<<")"<< endl;
  exit(1);

}



void pentagon(Topology* & Etopo, const double s12, const double s13, 
	      const double s14, const double s15, const double s23,
	      const double s24, const double s25, const double s34,
	      const double s35, const double s45, const double m1q,
	      const double m2q, const double m3q, const double m4q, 
	      const double m5q) {

  Matrix S(5,5);

  S.setValue(1,2,s12);
  S.setValue(2,1,s12);
  S.setValue(1,3,s13);
  S.setValue(3,1,s13);
  S.setValue(1,4,s14);
  S.setValue(4,1,s14);
  S.setValue(1,5,s15);
  S.setValue(5,1,s15);

  S.setValue(2,3,s23);
  S.setValue(3,2,s23);
  S.setValue(2,4,s24);
  S.setValue(4,2,s24);
  S.setValue(2,5,s25);
  S.setValue(5,2,s25);

  S.setValue(3,4,s34);
  S.setValue(4,3,s34);
  S.setValue(3,5,s35);
  S.setValue(5,3,s35);

  S.setValue(4,5,s45);
  S.setValue(5,4,s45);

  S.setValue(1,1,-2.0*m1q);
  S.setValue(2,2,-2.0*m2q);
  S.setValue(3,3,-2.0*m3q);
  S.setValue(4,4,-2.0*m4q);
  S.setValue(5,5,-2.0*m5q);

  Etopo = Topology::factory(S);

} 
void hexagon(Topology* & Ftopo, 
	     const double s12, const double s13, const double s14, 
	     const double s15, const double s16, 
	     const double s23, const double s24, const double s25, 
	     const double s26,
	     const double s34, const double s35, const double s36, 
	     const double s45, const double s46,
	     const double s56,
	     const double m1q, const double m2q, const double m3q, 
	     const double m4q, const double m5q, const double m6q) {

  Matrix S(6,6);

  S.setValue(1,2,s12);
  S.setValue(2,1,s12);
  S.setValue(1,3,s13);
  S.setValue(3,1,s13);
  S.setValue(1,4,s14);
  S.setValue(4,1,s14);
  S.setValue(1,5,s15);
  S.setValue(5,1,s15);
  S.setValue(1,6,s16);
  S.setValue(6,1,s16);

  S.setValue(2,3,s23);
  S.setValue(3,2,s23);
  S.setValue(2,4,s24);
  S.setValue(4,2,s24);
  S.setValue(2,5,s25);
  S.setValue(5,2,s25);
  S.setValue(2,6,s26);
  S.setValue(6,2,s26);

  S.setValue(3,4,s34);
  S.setValue(4,3,s34);
  S.setValue(3,5,s35);
  S.setValue(5,3,s35);
  S.setValue(3,6,s36);
  S.setValue(6,3,s36);

  S.setValue(4,5,s45);
  S.setValue(5,4,s45);
  S.setValue(4,6,s46);
  S.setValue(6,4,s46);

  S.setValue(5,6,s56);
  S.setValue(6,5,s56);


  S.setValue(1,1,-2.0*m1q);
  S.setValue(2,2,-2.0*m2q);
  S.setValue(3,3,-2.0*m3q);
  S.setValue(4,4,-2.0*m4q);
  S.setValue(5,5,-2.0*m5q);
  S.setValue(6,6,-2.0*m6q);

  Ftopo = Topology::factory(S);

} 

double Topology::getZ(  int i) const {

  assert(( i <= mydata->getNpoint() ) && ( (mydata->z) != 0 )); 
  return ((mydata->z)[i-1]);
}

double Topology::getRange( int i) const {

  assert( ( i <= mydata->getNpoint() ) && ( ( mydata->range) != 0 ));
  return mydata->range[i-1];
}

Topology::~Topology(){

  /*
   * Clear the vector with the references to the subtopologies,
   * and the map where we store evaluated integrals for this topology.
   * Note that we do not need to delete the contents. This
   * will be done when data::clear is called.
   */
  cache.clear();
  subtopos.clear();

}

Topology* Topology::factory(const Matrix& Sin){
 
  data *tmpdata = data::lookup(Sin);
  
  if ( tmpdata->getParent() != 0 ) {
    return( tmpdata->getParent() );
  }

  const int npoint = tmpdata->getNpoint(); 

  Topology* tmptopo = new Topology();
  tmptopo->mydata = tmpdata;
  tmpdata->setParent(tmptopo);

  if ( npoint  > 2 ) {
    for (int i=1; i <= npoint; i++){
      tmptopo->subtopos.push_back( tmptopo->evalSubtopo(i) );	
    }
  } 

  return(tmptopo);

}


void Topology::cacheclear(){

  data::clear();

}

Topology* Topology::evalSubtopo( int i) {

  return( factory( (*(mydata->S)).getSubMatrix(i) ) );

}



/*data dbox(FourMomentum q1, FourMomentum q2, FourMomentum q3,
		    FourMomentum q4, 
		    double m1, double m2, double m3, double m4){
  vector <FourMomentum> S1;
  vector <double> S1m;
  S1.push_back(q1);
  S1.push_back(q2);
  S1.push_back(q3);
  S1.push_back(q4);
  S1m.push_back(m1);  
  S1m.push_back(m2);
  S1m.push_back(m3);
  S1m.push_back(m4);
  data tmp(S1,S1m);
  return tmp;
  }*/

