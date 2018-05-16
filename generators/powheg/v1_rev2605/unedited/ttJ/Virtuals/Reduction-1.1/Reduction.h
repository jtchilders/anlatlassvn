/* $Modified: Wed Feb  3 09:36:50 2010 by uwer $ */
#ifndef _REDUCTION_H_
#define _REDUCTION_H_
#include "FourMomentum.h"
#include <vector>
#include <cstring>
#include <iostream>
#include <complex>



#define ONLYREALPART

#ifdef ONLYREALPART
typedef double IntType;
inline IntType mycast(std::complex<double> cint) {
  return cint.real();
}

#else
typedef std::complex<double> IntType;
inline IntType mycast(std::complex<double> cint) {
  return cint;
}
#endif







const double tiny2 = 1.e-7;

enum Bs {B0=0,B1=1,B00=2,B11=3,MAXB=4};

enum XBs {XB0=B0,XB1=B1,XB22=B00,XB21=B11};

enum Cs {C0=0,C1=1,C2=2,
	 C00=3,C11=4,C22=5,C12=6,
	 C001=7,C002=8,C111=9,C222=10,C112=11,C122=12,MAXC=13};

enum XCs {XC0=C0, XC11=C1,XC12=C2,
	  XC24=C00,XC21=C11,XC22=C22,XC23=C12,
	  XC35=C001,XC36=C002,XC31=C111,XC32=C222,XC33=C112,XC34=C122};


enum Ds {D0=0,D1=1,D2=2,D3=3,
	 D00=4,D11=5,D22=6,D33=7,D12=8,D13=9,D23=10,
	 D001=11,D002=12,D003=13,D111=14,D222=15,D333=16,D112=17,
	 D113=18,D122=19,D223=20,D133=21,D233=22,D123=23,
	 D0000=24,D0011=25,D0022=26,D0033=27,D0012=28,D0013=29,D0023=30,
	 D1111=31,D2222=32,D3333=33,D1112=34,D1113=35,D1222=36,D2223=37,
	 D1333=38,D2333=39,D1122=40,D1133=41,D2233=42,D1123=43,
	 D1223=44,D1233=45,MAXD=46};

enum XDs {XD0=D0,XD11=D1,XD12=D2,XD13=D3,
	  XD27=D00,XD21=D11,XD22=D22,XD23=D33,XD24=D12,XD25=D13,XD26=D23,
	  XD311=D001,XD312=D002,XD313=D003,XD31=D111,XD32=D222,XD33=D333,
	  XD34=D112,
	  XD35=D113,XD36=D122,XD38=D223,XD37=D133,XD39=D233,XD310=D123,
	  XD422=D0000,XD416=D0011,XD417=D0022,XD418=D0033,
	  XD419=D0012,XD420=D0013,XD421=D0023,
	  XD41=D1111,XD42=D2222,XD43=D3333,XD44=D1112,XD45=D1113,XD46=D1222,
	  XD47=D2223,
	  XD48=D1333,XD49=D2333,XD410=D1122,XD411=D1133,XD412=D2233,XD413=D1123,
	  XD414=D1223,XD415=D1233};


class Dinput {
 private:
  int masscount;
  int external_masscount;
  void countmasses(double m);
  void countextmasses(double mq);

  //  FourMomentum p1,p2,p3,p5,p6,p7;

 public:
  double p1q,p2q,p3q,p5q,p6q,p7q,p1p2,p1p3,p2p3;
  double m0q,m1q,m2q,m3q;
  Dinput(const double & k1q, const double & k2q, 
  	 const double & k3q, const double & k4q,
  	 const double & k5q, const double & k7q,
  	 const double tm0q, const double tm1q, const double tm2q, 
  	 const double tm3q);

  void printmaple() const; 
  int getMassCount() const; 
  int getExtMassCount() const;
  //  ~Dinput(){}

};

class Cinput {

 private:
  int masscount;
  int external_masscount;

  void countmasses(double m);
  void countextmasses(double mq);
  //  FourMomentum p1,p2,p5;

 public:
  double p1q,p2q,p5q,p1p2;
  double m0q,m1q,m2q;
  //  Cinput(const FourMomentum & k1,const FourMomentum & k2,
  //	 const double tm0q, const double tm1q, const double tm2q){
  Cinput(const double k1q, const double k2q, const double k5q,
  	 const double tm0q, const double tm1q, const double tm2q);

  void printmaple() const;

  int getMassCount() const;
  int getExtMassCount() const; 
  //  ~Cinput(){}

};

class Binput {
 public:
  double m0q,m1q;
  double p1q;
  Binput(double k1q, const double tm0q, const double tm1q);
  void printmaple();
  //  ~Binput(){}
};

template <class Input> class Coeffs {
 private:	
  Input *para;
  IntType *res_ptr;	
 public:
  Coeffs(Input *tmppara, IntType *tmp_res_ptr){
    para = tmppara;
    res_ptr = tmp_res_ptr;
  }
  ~Coeffs(){
    delete para; 
    delete[] res_ptr;
  }
  Input *getpara() { return para; }
  IntType *getcoeffptr() { return res_ptr; }	
};


class CoeffCache {
  private:
  std::vector<Coeffs<Dinput>*> dcache;
  std::vector<Coeffs<Cinput>*> ccache;
  std::vector<Coeffs<Binput>*> bcache;

  public:
  ~CoeffCache(){
    reset();
  }
  void reset();
  
  void lookup(IntType* &ptr, 
	      const double & k1q, const double & k2q, 
	      const double & k3q, const double & k4q, 
	      const double & k5q, const double & k7q, 
	      const double m0q, const double m1q, const double m2q, 
	      const double m3q);

  void lookup(IntType* &ptr, const double k1q, const double k2q,
	      const double k5q,
	      const double m0q, const double m1q, const double m2q);

  void lookup(IntType* &ptr, const double k1q,
	      const double m0q, const double m1q);

  IntType* evalCoeff(const Dinput & in);
  IntType* evalCoeff(const Cinput & in);
  IntType* evalCoeff(const Binput & in);
};

#endif //_REDUCTION_H_
