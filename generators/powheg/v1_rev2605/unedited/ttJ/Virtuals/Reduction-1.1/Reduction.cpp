// $Modified: Fri Aug 11 10:07:54 2006 by puwer $
#include "Reduction.h"
#include "ScalarInt.h"
#include <cmath>
#include <algorithm>
#include <limits>

//#define REDUCTION_DEBUG

#define DEBUGPRINT(a_) cout << #a_ << " = "<< a_ << endl
#define PRINT(a_,b_) cout << #a_ << b_ << ";" << endl;

using namespace std;

const double expflag = 1.0;
 
const double tinymass = 1e-08;

Binput::Binput(double k1q, const double tm0q, const double tm1q){
  p1q = k1q;
  m0q = tm0q;
  m1q = tm1q;
}

void Binput::printmaple(){
  PRINT(p1q:=,p1q); 
  PRINT(m1:=,sqrt(m0q));
  PRINT(m2:=,sqrt(m1q));
}

void Cinput::countmasses(double m){
  if (fabs(m) > 1.0){ 
    masscount++;
  } else {
    if (fabs(m) != 0.0)
      std::cout << "Mass value " << fabs(m) << " encountered.";
  }
}

void Cinput::countextmasses(double mq){
  if (fabs(mq) > tiny2)
    external_masscount++;
}

Cinput::Cinput(const double k1q, const double k2q, const double k5q,
	       const double tm0q, const double tm1q, const double tm2q){
  //    p1 = k1;
  //p2 = k2;
  //p5 = p1 + p2;
  m0q = tm0q;
  m1q = tm1q;
  m2q = tm2q;

  /*
   * This is really odd....
   */
  masscount = 0;
  countmasses(m0q);
  countmasses(m1q);
  countmasses(m2q);

  external_masscount = 0;
  p1q = k1q;
  p2q = k2q;
  p5q = k5q;
  p1p2 = (k5q-k1q-k2q)/2.0;
  countextmasses(p1q);
  countextmasses(p2q);
  countextmasses(p5q);
}

void Cinput::printmaple() const {
  PRINT(p1q:=,p1q); 
  PRINT(p2q:=,p2q); 
  PRINT(p3q:=,p5q); 
  PRINT(p1p2:=,p1p2); 
  PRINT(m1:=,sqrt(m0q));
  PRINT(m2:=,sqrt(m1q));
  PRINT(m3:=,sqrt(m2q));
}

int Cinput::getMassCount() const {
  return(masscount);
} 

int Cinput::getExtMassCount() const {
  return(external_masscount);
}


void Dinput::countmasses(double m){
  if (m != 0.0) 
    masscount++;
}

void Dinput::countextmasses(double mq){
  
  if (fabs(mq) > tiny2) 
    external_masscount++;
  
}

Dinput::Dinput(const double & k1q, const double & k2q, 
	       const double & k3q, const double & k4q,
	       const double & k5q, const double & k7q,
	       const double tm0q, const double tm1q, const double tm2q, 
	       const double tm3q){
  //  Dinput(const FourMomentum & k1, const FourMomentum & k2, 
  //	 const FourMomentum & k3,
  //	 const double tm0q, const double tm1q, const double tm2q, 
  //	 const double tm3q){
    //   p1 = k1;
    //p2 = k2;
    //p3 = k3;
    //p5 = p1+p2;
    //p6 = p1+p2+p3;
    //p7 = p2+p3;

  p1q = k1q;
  p2q = k2q;
  p3q = k3q;
  p5q = k5q;
  p6q = k4q;
  p7q = k7q;
  p1p2 = (p5q - p1q - p2q) / 2.0;
  p1p3 = (p6q - p5q + p2q - p7q)/2.0;
  p2p3 = (p7q - p2q - p3q)/2.0;
  
  m0q = tm0q;
  m1q = tm1q;
  m2q = tm2q;
  m3q = tm3q;
  
  /*
   * This is really odd....
   */
  masscount = 0;
  countmasses(m0q);
  countmasses(m1q);
  countmasses(m2q);
  countmasses(m3q);
  
  external_masscount = 0;
  countextmasses(p1q);
  countextmasses(p2q);
  countextmasses(p3q);
  countextmasses(p5q);
  countextmasses(p6q);
  countextmasses(p7q);
}

void Dinput::printmaple() const {
    PRINT(p1q:=,p1q); 
    PRINT(p2q:=,p2q); 
    PRINT(p3q:=,p3q); 
    PRINT(p4q:=,p6q); 
    PRINT(p12q:=,p5q); 
    PRINT(p23q:=,p7q);  
    PRINT(p1p2:=,p1p2); 
    PRINT(p1p3:=,p1p3); 
    PRINT(p2p3:=,p2p3);
    PRINT(m0:=,sqrt(m0q));
    PRINT(m1:=,sqrt(m1q));
    PRINT(m2:=,sqrt(m2q));
    PRINT(m3:=,sqrt(m3q));
}

int Dinput::getMassCount() const {
  return(masscount);
} 

int Dinput::getExtMassCount() const {
  return(external_masscount);
}



void CoeffCache::reset(){
  //  for (unsigned int i=0; i < dcache.size(); i++) delete dcache[i];
  for ( vector<Coeffs<Dinput>*>::iterator i = dcache.begin();
	i != dcache.end(); i++ ) {
    delete *i;
  }
  dcache.clear();

  //for (unsigned int i=0; i < ccache.size(); i++) delete ccache[i];
  for ( vector<Coeffs<Cinput>*>::iterator i = ccache.begin();
	i != ccache.end(); i++ ) {
    delete *i;
  }
  ccache.clear();

  // for (unsigned int i=0; i < bcache.size(); i++) delete bcache[i];
  for ( vector<Coeffs<Binput>*>::iterator i = bcache.begin();
	i != bcache.end(); i++ ) {
    delete *i;
  }
  bcache.clear();
  //   std::cout << "cache cleared..." << std::endl; 
}


void CoeffCache::lookup(IntType* &ptr, 
			const double & k1q, const double & k2q, 
			const double & k3q, const double & k4q, 
			const double & k5q, const double & k7q, 
			const double m0q, const double m1q, const double m2q, 
			const double m3q){
    
  Dinput* tmpinput=new Dinput(k1q,k2q,k3q,k4q,k5q,k7q,m0q,m1q,m2q,m3q);
  
  //  for( unsigned int i=0; i< dcache.size();i++){
  for ( vector<Coeffs<Dinput>*>::iterator i = dcache.begin();
	i != dcache.end(); i++ ) {
    bool exist = 
      (0 == memcmp( (*i)->getpara(),tmpinput,sizeof(Dinput)));
    if (exist) {
      // std::cout << "fund in dcache..." << std::endl;
      delete tmpinput;
      ptr = (*i)->getcoeffptr();
	return;
    }
  }

  // std::cout << "not fund, calculate it..." << std::endl;

  Coeffs<Dinput> *newcoeff = 
    new Coeffs<Dinput>(tmpinput, evalCoeff(*tmpinput));

  dcache.push_back(newcoeff);

  ptr = (*newcoeff).getcoeffptr();

}

void CoeffCache::lookup(IntType* &ptr, const double k1q, const double k2q,
			const double k5q,
			const double m0q, const double m1q, const double m2q){
  
  Cinput* tmpinput = new Cinput(k1q,k2q,k5q,m0q,m1q,m2q);
  
  //  for( unsigned int i=0; i< ccache.size();i++){
  for ( vector<Coeffs<Cinput>*>::iterator i = ccache.begin();
	i != ccache.end(); i++ ) {
    bool exist = 
      (0==memcmp( (*i)->getpara(),tmpinput,sizeof(Cinput)));
    if (exist) {
      //	std::cout << "fund in dcache..." << std::endl;
      delete tmpinput;
      ptr = (*i)->getcoeffptr();
      return;
    }
  }

  //    std::cout << "not fund, calculate it..." << std::endl;

  Coeffs<Cinput> *newcoeff = 
    new Coeffs<Cinput>(tmpinput, evalCoeff(*tmpinput));

  ccache.push_back(newcoeff);

  ptr = (*newcoeff).getcoeffptr();

}

void CoeffCache::lookup(IntType* &ptr, const double k1q,
			const double m0q, const double m1q){

  Binput* tmpinput = new Binput(k1q,m0q,m1q);

  
  //  for( unsigned int i=0; i< bcache.size();i++){
  for ( vector<Coeffs<Binput>*>::iterator i = bcache.begin();
	i != bcache.end(); i++ ) {
    bool exist = 
      (0==memcmp( (*i)->getpara(),tmpinput,sizeof(Binput)));
    if (exist) {
      //	std::cout << "fund in dcache..." << std::endl;
      delete tmpinput;
      ptr = (*i)->getcoeffptr();
      return;
    }
  }

  //    std::cout << "not fund, calculate it..." << std::endl;

  Coeffs<Binput> *newcoeff = 
    new Coeffs<Binput>(tmpinput, evalCoeff(*tmpinput));

  bcache.push_back(newcoeff);

  ptr = (*newcoeff).getcoeffptr();

}

IntType*  CoeffCache::evalCoeff(const Binput & in){
  IntType *ptr = new IntType[MAXB];
  double p1q = in.p1q;

  //----------------------------  
  // fill(ptr,ptr+MAXB,1.0);
  // return(ptr);
  //----------------------------  

  if ( abs(p1q) > tinymass ) {

    ptr[B0] = ScalarInt(in);

    ptr[XB1] = 0.5/p1q * ( ScalarInt(in.m0q) - ScalarInt(in.m1q) 
			   + ( in.m1q - in.m0q - p1q ) * ptr[B0] );

    ptr[XB21] = 1.0 / 6.0 / p1q * 
      ( 2.0 * ScalarInt(in.m1q) - 2.0 * in.m0q * ptr[B0] 
	- 4.0 * (in.m0q - in.m1q + p1q) * ptr[XB1]  
	- ( in.m0q + in.m1q - 1.0/3.0*p1q ) * expflag );
 
    ptr[XB22] = 1.0/6.0 * 
      ( ScalarInt(in.m1q) + 2.0 * in.m0q*ptr[B0]
	+ (in.m0q-in.m1q+p1q ) * ptr[XB1] 
	+ ( in.m0q + in.m1q - 1.0/3.0*p1q ) * expflag );   

    return ptr;
  };

  if ((abs(in.m0q)<tinymass) && (abs(in.m1q)<tinymass) ) {
    /* scale less integral set everything to zero
       ptr[B0] = 0.0;
       ptr[XB1]  = 0.0;
       ptr[XB22] = 0.0;
       ptr[XB21] = 0.0;
    */
    ptr[B0] = ScalarInt(in);
    ptr[XB1]  = -1.0/2.0 *  ptr[B0];
    //ptr[XB22] = 0.0; // ori 
    ptr[XB22] = 1.0/4.0 * ScalarInt(0.0); //changed 06.06.06
    ptr[XB21] = 1.0/3.0 * ptr[B0];
    return ptr;
  };

  /*
   *  Check if the masses are equal:
   */
  
  if ( abs(in.m0q/in.m1q - 1.0 ) < 0.00001 ) {
    ptr[B0] = ScalarInt(in);
    ptr[XB1] = -1.0/2.0 * ScalarInt(in.m0q) / in.m0q + 1.0/2.0 * expflag ;
    ptr[XB22] = 1.0/2.0 * ScalarInt(in.m0q) ;
    ptr[XB21] = 1.0/3.0 * ScalarInt(in.m0q) / in.m0q - 1.0/3.0 * expflag ; 
    return ptr;

  } else {
    /* 
     * pq = 0, m0 \not = m1 
     */
    cout << "not yet checked...\n";
    exit(1);
    //ptr[B0] = (ScalarInt(in.m0q)-ScalarInt(in.m1q))/(in.m0q-in.m1q);
    ptr[B0] = ScalarInt(in);
    ptr[XB21] = -1.0/3.0/pow(in.m0q-in.m1q,3) * 
      ( ( 3.0 * in.m0q*in.m0q + in.m1q*in.m1q - 3.0 * in.m0q*in.m1q) 
	* ScalarInt(in.m1q) - in.m0q*in.m0q * ScalarInt(in.m0q) )
      + (5.0 * in.m0q*in.m0q + 5.0*in.m0q*in.m1q - 4.0 * in.m1q*in.m1q)/18.0
      /(in.m0q-in.m1q)/(in.m0q-in.m1q) * expflag;
    
    ptr[XB1] = -1.0/4.0/pow(in.m0q-in.m1q,2) *
      ( (2.0*in.m1q-4.0*in.m0q)*ScalarInt(in.m1q) 
	+ 2.0 * in.m0q * ScalarInt(in.m0q) )
      - 1.0/4.0 * (in.m0q+in.m1q)/(in.m0q-in.m1q) * expflag;
    
    ptr[XB22] = 1.0/4.0 * ( 1.0/2.0* ( in.m1q*ScalarInt(in.m1q)
				       -in.m0q*ScalarInt(in.m0q))
			    / (in.m1q-in.m0q)
			    - (in.m1q-in.m0q)*ptr[XB1] + in.m0q * ptr[B0])
      + 3.0/16.0 * ( in.m0q + in.m1q ) * expflag;
    
    return ptr;
  };
 
  // We should never reach this point
  cout << "CoeffCache::evalCoeff(const Binput & ): Alien integral detected" 
       << endl << "Don't know what to do!" << endl;

  exit(1);

}

static void solve2(double Ainv[2][2], IntType ptr[],int a0,int a1, 
		     const IntType &  b0, const IntType & b1){
  ptr[a0] = Ainv[0][0] * b0 + Ainv[0][1] * b1;
  ptr[a1] = Ainv[1][0] * b0 + Ainv[1][1] * b1;
}

static void solve3(double Ainv[3][3], IntType ptr[],int a0,int a1,int a2, 
		   const IntType &  b0, const IntType & b1, const IntType & b2)
{

  ptr[a0] = Ainv[0][0] * b0 + Ainv[0][1] * b1 + Ainv[0][2] * b2;
  ptr[a1] = Ainv[1][0] * b0 + Ainv[1][1] * b1 + Ainv[1][2] * b2;
  ptr[a2] = Ainv[2][0] * b0 + Ainv[2][1] * b1 + Ainv[2][2] * b2;
}


IntType*  CoeffCache::evalCoeff(const Cinput & in){
  IntType *ptr = new IntType[MAXC];
  IntType *b12_ptr, *b13_ptr, *b23_ptr;
  IntType R1,R2,R3,R4,R5,R6,R8,R9,R10,R11,R12,R13,R14,R15;
  double f1,f2; 
  // double p1q = dotp(in.p1,in.p1); 
  // double p1p2 = dotp(in.p1,in.p2);
  // double p2q = dotp(in.p2,in.p2);
  double p1q  = in.p1q; 
  double p2q  = in.p2q;
  double p1p2 = in.p1p2;
  double p5q = in.p5q;

  double X2inv[2][2],det2;

  //----------------------------  
  // fill(ptr,ptr+MAXC,1.0);
  // return(ptr);
  //----------------------------  

  // When setting the scalar-integrals to very large values
  // massless momenta which are not exactly massless create
  // finite contributions. To avoid this we set the squared
  // momenta to zero...
  if (fabs(p1q)<tinymass)
    p1q = 0.0;

  if (fabs(p2q)<tinymass)
    p2q = 0.0;

  if (fabs(p1p2)<tinymass)
    p1p2 = 0.0;


  det2 = p1q*p2q-p1p2*p1p2;

#define CHECKINVGRAMDET
#ifdef CHECKINVGRAMDET
  if (det2 == 0.0 ) {
    cout<<"Problem inverse Gram det=0 in CoeffCache::evalCoeff(const Cinput & in)"<<endl;
    cout.flush();
    det2=numeric_limits<double>::quiet_NaN();
  }
#endif
  X2inv[0][0] = p2q/det2;
  X2inv[0][1] = X2inv[1][0] = -p1p2/det2;
  X2inv[1][1] = p1q/det2;


  f1 = in.m1q - in.m0q - p1q;
  f2 = in.m2q - in.m1q - 2.0 * p1p2 - p2q;
  
  lookup(b12_ptr,p1q,in.m0q,in.m1q);
  lookup(b13_ptr,p5q,in.m0q,in.m2q);
  lookup(b23_ptr,p2q,in.m1q,in.m2q);

  ptr[C0] = ScalarInt(in);

  R1 = 0.5 * ( f1*ptr[C0] + b13_ptr[B0] - b23_ptr[B0] );
  R2 = 0.5 * ( f2*ptr[C0] + b12_ptr[B0] - b13_ptr[B0] );

  solve2(X2inv,ptr,XC11,XC12,R1,R2);


  ptr[XC24] = 0.5 * (in.m0q*ptr[C0] 
		     + 0.5 * ( b23_ptr[B0] - f1 * ptr[XC11] 
			       - f2 * ptr[XC12] ) ) + 1.0/4.0 * expflag ;

  R3 = 0.5 * (f1*ptr[XC11] + b13_ptr[XB1] + b23_ptr[B0]) - ptr[XC24];
  R4 = 0.5 * (f1*ptr[XC12] + b13_ptr[XB1] - b23_ptr[XB1]);
  R5 = 0.5 * (f2*ptr[XC11] + b12_ptr[XB1] - b13_ptr[XB1]);
  R6 = 0.5 * (f2*ptr[XC12] - b13_ptr[XB1] ) - ptr[XC24];

  solve2(X2inv,ptr,XC21,XC23,R3,R5);

#ifdef REDUCTION_DEBUG
  IntType tmpC23 = ptr[XC23];
#endif

  solve2(X2inv,ptr,XC23,XC22,R4,R6);

#ifdef REDUCTION_DEBUG
  cout.precision(15);
  if ( abs(ptr[XC23]/tmpC23-1.0) > 0.000001 ) 
    cout << "Got two different values for C23:"
	 << ptr[XC23] << "," << tmpC23 << endl;
#endif

  R11 = 0.5 * ( f1 * ptr[XC24] + b13_ptr[XB22] - b23_ptr[XB22] );
  R15 = 0.5 * ( f2 * ptr[XC24] + b12_ptr[XB22] - b13_ptr[XB22]) ;

  solve2(X2inv,ptr,XC35,XC36,R11,R15);


  R8  = 0.5 * ( f1 * ptr[XC21] + b13_ptr[XB21] - b23_ptr[B0] ) 
    - 2.0 * ptr[XC35];
  R9  = 0.5 * ( f1 * ptr[XC22] + b13_ptr[XB21] - b23_ptr[XB21] );
  R10 = 0.5 * ( f1 * ptr[XC23] + b13_ptr[XB21] + b23_ptr[XB1] ) - ptr[XC36];

  R12 = 0.5 * ( f2 * ptr[XC21] + b12_ptr[XB21] - b13_ptr[XB21] );
  R13 = 0.5 * ( f2 * ptr[XC22] - b13_ptr[XB21] ) - 2.0 * ptr[XC36];
  R14 = 0.5 * ( f2 * ptr[XC23] - b13_ptr[XB21] ) - ptr[XC35];

  solve2(X2inv,ptr,XC31,XC33,R8,R12);
  solve2(X2inv,ptr,XC34,XC32,R9,R13);

#ifdef REDUCTION_DEBUG
  IntType tmpC33 = ptr[XC33];
  IntType tmpC34 = ptr[XC34];
#endif

  solve2(X2inv,ptr,XC33,XC34,R10,R14);

#ifdef REDUCTION_DEBUG
  if ( abs(ptr[XC33]/tmpC33-1.0) > 0.000001 ) 
    cout << "Got two different values for C33:"
	 << ptr[XC33] << "," << tmpC33 << endl;
  if ( abs(ptr[XC34]/tmpC34-1.0) > 0.000001 ) 
    cout << "Got two different values for C34:"
	      << ptr[XC34] << "," << tmpC34 << endl;
#endif

  return ptr;

}

IntType*  CoeffCache::evalCoeff(const Dinput & in){
  IntType *ptr = new IntType[MAXD];
  IntType *c123_ptr, *c124_ptr,*c134_ptr,*c234_ptr;
  IntType R20,R21,R22,R30,R31,R32,R33,R34,R35,R36,R37,R38;
  IntType R41,R42,R43,R44,R45,R46,R50,R51,R52,R53,R54,R55,R56,R57,R58;
  IntType R60,R61,R62,R63,R64,R65,R66,R67,R68,R69,R70,R71,R72,R73,R74,R75,
    R76,R77;
  double Xinv[3][3];
  // double p1q = dotp(in.p1,in.p1);
  // double p2q = dotp(in.p2,in.p2);
  // double p3q = dotp(in.p3,in.p3);
  // double p5q = dotp(in.p5,in.p5);
  // double p6q = dotp(in.p6,in.p6);
  // double p1p2 = dotp(in.p1,in.p2);
  // double p1p3 = dotp(in.p1,in.p3);
  // double p2p3 = dotp(in.p2,in.p3);  

  const double p1q = in.p1q;
  const double p2q = in.p2q;
  const double p3q = in.p3q;
  const double p5q = in.p5q;
  const double p6q = in.p6q;
  const double p7q = in.p7q;
  const double p1p2 = in.p1p2;
  const double p1p3 = in.p1p3;
  const double p2p3 = in.p2p3;
  
  double det3;
  double f1,f2,f3;
  

  //----------------------------  
  // fill(ptr,ptr+MAXD,1.0);
  // return(ptr);
  //----------------------------  

  //  lookup(c123_ptr,in.p1,in.p2,in.m0q,in.m1q,in.m2q);
  //  lookup(c124_ptr,in.p1,in.p7,in.m0q,in.m1q,in.m3q);
  //  lookup(c134_ptr,in.p5,in.p3,in.m0q,in.m2q,in.m3q);
  //  lookup(c234_ptr,in.p2,in.p3,in.m1q,in.m2q,in.m3q);

  lookup(c123_ptr,p1q,p2q,p5q,in.m0q,in.m1q,in.m2q);
  lookup(c124_ptr,p1q,p7q,p6q,in.m0q,in.m1q,in.m3q);
  lookup(c134_ptr,p5q,p3q,p6q,in.m0q,in.m2q,in.m3q);
  lookup(c234_ptr,p2q,p3q,p7q,in.m1q,in.m2q,in.m3q);

  f1 = in.m1q - in.m0q - p1q;
  f2 = in.m2q - in.m1q + p1q - p5q;
  f3 = in.m3q - in.m2q + p5q - p6q;

  Xinv[0][0] = p2q*p3q - p2p3*p2p3;
  Xinv[0][1] = - ( p1p2*p3q-p1p3*p2p3 );
  Xinv[0][2] = ( p1p2*p2p3-p1p3*p2q );

  det3 = p1q*Xinv[0][0] + p1p2*Xinv[0][1] + p1p3 * Xinv[0][2];
  
#ifdef REDUCTION_DEBUG
  DEBUGPRINT(det3);
#endif

  Xinv[0][0] = Xinv[0][0]/det3;
  Xinv[0][1] = Xinv[0][1]/det3;
  Xinv[0][2] = Xinv[0][2]/det3;

  Xinv[1][1] = ( p1q*p3q-p1p3*p1p3)/det3;
  Xinv[1][2] = - ( p1q*p2p3-p1p3*p1p2)/det3; 
  Xinv[2][2] = (p1q*p2q-p1p2*p1p2)/det3;

  Xinv[1][0] = Xinv[0][1];
  Xinv[2][0] = Xinv[0][2];
  Xinv[2][1] = Xinv[1][2];

#ifdef REDUCTION_DEBUG
  DEBUGPRINT(Xinv[0][0]);
  DEBUGPRINT(Xinv[0][1]);
  DEBUGPRINT(Xinv[1][0]);
  DEBUGPRINT(Xinv[0][2]);
  DEBUGPRINT(Xinv[2][0]);
  DEBUGPRINT(Xinv[1][1]);
  DEBUGPRINT(Xinv[1][2]);
  DEBUGPRINT(Xinv[2][1]);
  DEBUGPRINT(Xinv[2][2]);
#endif

  ptr[D0] = ScalarInt(in);

  R20 = 0.5 * ( f1 * ptr[D0] + c134_ptr[C0] - c234_ptr[C0] );
  R21 = 0.5 * ( f2 * ptr[D0] + c124_ptr[C0] - c134_ptr[C0] );
  R22 = 0.5 * ( f3 * ptr[D0] + c123_ptr[C0] - c124_ptr[C0] );

  solve3(Xinv,ptr,XD11,XD12,XD13,R20,R21,R22);

  ptr[XD27] = in.m0q*ptr[D0]-0.5*( f1 * ptr[XD11] + f2 * ptr[XD12]
				   + f3 * ptr[XD13] - c234_ptr[C0]);

  R30 = 0.5 * ( f1 * ptr[XD11] + c134_ptr[XC11] + c234_ptr[C0] ) - ptr[XD27];
  R31 = 0.5 * ( f2 * ptr[XD11] + c124_ptr[XC11] - c134_ptr[XC11] );
  R32 = 0.5 * ( f3 * ptr[XD11] + c123_ptr[XC11] - c124_ptr[XC11] );

  R33 = 0.5 * ( f1 * ptr[XD12] + c134_ptr[XC11] - c234_ptr[XC11]); 
  R34 = 0.5 * ( f2 * ptr[XD12] + c124_ptr[XC12] - c134_ptr[XC11]) - ptr[XD27]; 
  R35 = 0.5 * ( f3 * ptr[XD12] + c123_ptr[XC12] - c124_ptr[XC12]); 

  R36 = 0.5 * ( f1 * ptr[XD13] + c134_ptr[XC12] - c234_ptr[XC12]); 
  R37 = 0.5 * ( f2 * ptr[XD13] + c124_ptr[XC12] - c134_ptr[XC12]); 
  
  R38 = 0.5 * ( f3 * ptr[XD13] - c124_ptr[XC12] ) - ptr[XD27]; 

  solve3(Xinv,ptr,XD21,XD24,XD25,R30,R31,R32);
  solve3(Xinv,ptr,XD24,XD22,XD26,R33,R34,R35);
  solve3(Xinv,ptr,XD25,XD26,XD23,R36,R37,R38);

  ptr[XD311] = 0.5*(in.m0q*ptr[XD11]
		    - 0.5 * ( f1 * ptr[XD21] + f2 * ptr[XD24] 
			      + f3 * ptr[XD25] + c234_ptr[C0]));

  ptr[XD312] = 0.5*(in.m0q*ptr[XD12]
		    - 0.5 * ( f1 * ptr[XD24] + f2 * ptr[XD22] 
			      + f3 * ptr[XD26] - c234_ptr[XC11]));

  ptr[XD313] = 0.5*(in.m0q*ptr[XD13]
		    - 0.5 * ( f1 * ptr[XD25] + f2 * ptr[XD26] 
			      + f3 * ptr[XD23] - c234_ptr[XC12]));

  R41 = 0.5 * ( f1 * ptr[XD21] - c234_ptr[C0] + c134_ptr[XC21]) 
    - 2.0 * ptr[XD311];
  R42 = 0.5 * ( f2 * ptr[XD21] - c134_ptr[XC21] + c124_ptr[XC21]); 
  R43 = 0.5 * ( f3 * ptr[XD21] - c124_ptr[XC21] + c123_ptr[XC21]); 

  R44 = 0.5 * ( f1 * ptr[XD24] + c234_ptr[XC11] + c134_ptr[XC21]) 
    - ptr[XD312];
  R45 = 0.5 * ( f2 * ptr[XD24] - c134_ptr[XC21] + c124_ptr[XC23])
    - ptr[XD311];
  R46 = 0.5 * ( f3 * ptr[XD24] - c124_ptr[XC23] + c123_ptr[XC23]);

  R50 = 0.5 * ( f1 * ptr[XD22] - c234_ptr[XC21] + c134_ptr[XC21]);

  R51 = 0.5 * ( f2 * ptr[XD22] - c134_ptr[XC21] + c124_ptr[XC22])
    -2.0 * ptr[XD312];

  R52 = 0.5 * ( f3 * ptr[XD22] - c124_ptr[XC22] + c123_ptr[XC22]);

  R56 = 0.5 * ( f1 * ptr[XD23] - c234_ptr[XC22] + c134_ptr[XC22]);

  R57 = 0.5 * ( f2 * ptr[XD23] - c134_ptr[XC22] + c124_ptr[XC22]);

  R58 = 0.5 * ( f3 * ptr[XD23] - c124_ptr[XC22]) - 2.0 * ptr[XD313];

  solve3(Xinv,ptr,XD31,XD34,XD35,R41,R42,R43);
  solve3(Xinv,ptr,XD36,XD32,XD38,R50,R51,R52);
  solve3(Xinv,ptr,XD37,XD39,XD33,R56,R57,R58);
  solve3(Xinv,ptr,XD34,XD36,XD310,R44,R45,R46);

  ptr[XD416] = 1.0/3.0 * ( in.m0q * ptr[XD21] 
			   + 0.5 * (c234_ptr[C0] - f1 * ptr[XD31] 
				    - f2 * ptr[XD34] - f3 * ptr[XD35] ) );

  ptr[XD417] = 1.0/3.0 * ( in.m0q * ptr[XD22] 
			   + 0.5 * (c234_ptr[XC21] - f1 * ptr[XD36] 
				    - f2 * ptr[XD32] - f3 * ptr[XD38] ) );

  ptr[XD418] = 1.0/3.0 * ( in.m0q * ptr[XD23] 
			   + 0.5 * (c234_ptr[XC22] - f1 * ptr[XD37] 
				    - f2 * ptr[XD39] - f3 * ptr[XD33] ) );

  ptr[XD419] = 1.0/3.0 * ( in.m0q * ptr[XD24] 
			   - 0.5 * (c234_ptr[XC11] + f1 * ptr[XD34] 
				    + f2 * ptr[XD36] + f3 * ptr[XD310] ) );

  ptr[XD420] = 1.0/3.0 * ( in.m0q * ptr[XD25] 
			   - 0.5 * (c234_ptr[XC12] + f1 * ptr[XD35] 
				    + f2 * ptr[XD310] + f3 * ptr[XD37] ) );

  ptr[XD421] = 1.0/3.0 * ( in.m0q * ptr[XD26] 
			   + 0.5 * (c234_ptr[XC23] - f1 * ptr[XD310] 
				    - f2 * ptr[XD38] - f3 * ptr[XD39] ) );

  ptr[XD422] = 1.0/6.0 * ( in.m0q * ptr[XD27] + c234_ptr[XC24] 
			   - p1q * ptr[XD416] 
			   - p2q * ptr[XD417]
			   - p3q * ptr[XD418] 
			   - 2.0 * p1p2 * ptr[XD419]
			   - 2.0 * p1p3 * ptr[XD420] 
			   - 2.0 * p2p3 * ptr[XD421] + 1.0/12.0 * expflag);

  R60 = 0.5 * ( f1 * ptr[XD32] - c234_ptr[XC31] + c134_ptr[XC31] );
  R61 = 0.5 * ( f2 * ptr[XD32] - c134_ptr[XC31] + c124_ptr[XC32] ) 
    - 3.0*ptr[XD417];
  R62 = 0.5 * ( f3 * ptr[XD32] - c124_ptr[XC32] + c123_ptr[XC32] );

  R63 = 0.5 * ( f1 * ptr[XD31] + c234_ptr[C0] + c134_ptr[XC31] ) 
    - 3.0*ptr[XD416];
  R64 = 0.5 * ( f2 * ptr[XD31] - c134_ptr[XC31] + c124_ptr[XC31] ); 
  R65 = 0.5 * ( f3 * ptr[XD31] - c124_ptr[XC31] + c123_ptr[XC31] ); 

  R66 = 0.5 * ( f1 * ptr[XD33] - c234_ptr[XC32] + c134_ptr[XC32] ); 
  R67 = 0.5 * ( f2 * ptr[XD33] - c134_ptr[XC32] + c124_ptr[XC32] ); 
  R68 = 0.5 * ( f3 * ptr[XD33] - c124_ptr[XC32] ) - 3.0 * ptr[XD418]; 

  R69 = 0.5 * ( f1 * ptr[XD37] + c234_ptr[XC22] + c134_ptr[XC34] )-ptr[XD418]; 
  R70 = 0.5 * ( f2 * ptr[XD37] - c134_ptr[XC34] + c124_ptr[XC34] ); 
  R71 = 0.5 * ( f3 * ptr[XD37] - c124_ptr[XC34] ) - 2.0 * ptr[XD420]; 

  R72 = 0.5 * ( f1 * ptr[XD34] - c234_ptr[XC11] + c134_ptr[XC31] )
    -2.0*ptr[XD419]; 
  R73 = 0.5 * ( f2 * ptr[XD34] - c134_ptr[XC31] + c124_ptr[XC33] )
    - ptr[XD416]; 
  R74 = 0.5 * ( f3 * ptr[XD34] - c124_ptr[XC33] + c123_ptr[XC33] ); 

  R75 = 0.5 * ( f1 * ptr[XD38] - c234_ptr[XC33] + c134_ptr[XC33] ); 
  R76 = 0.5 * ( f2 * ptr[XD38] - c134_ptr[XC33] + c124_ptr[XC32] )
    - 2.0*ptr[XD421]; 
  R77 = 0.5 * ( f3 * ptr[XD38] - c124_ptr[XC32] ) - ptr[XD417]; 

  solve3(Xinv,ptr,XD46,XD42,XD47,R60,R61,R62);
  solve3(Xinv,ptr,XD41,XD44,XD45,R63,R64,R65);
  solve3(Xinv,ptr,XD48,XD49,XD43,R66,R67,R68);
  solve3(Xinv,ptr,XD411,XD415,XD48,R69,R70,R71);
  solve3(Xinv,ptr,XD44,XD410,XD413,R72,R73,R74);
  solve3(Xinv,ptr,XD414,XD47,XD412,R75,R76,R77);

  return ptr;

}
