// $Modified: Thu Sep 21 12:08:13 2006 by puwer $
#include "Reduction.h"
#include "Evalf.h"
#include "ScalarInt.h"
#include <complex>
#include <iostream>

using namespace std;

double den(double x);

double den(double x){
  return(1.0/x);
}


void RemainingDiags_gg(IntType f_nfl[12][133], 
		    const FourMomentum & p1, const FourMomentum & p2, 
		    const FourMomentum & p3,
		    const FourMomentum & kq, const FourMomentum & kqb,
		    IntType* Bptr[], IntType* Cptr[],
		    IntType* Dptr[],
		    double rterms
		    ) {

  double m = kq.getmass();
  double s = 2.0 * dotp(p1,p2);
  double y  = 2.0 * dotp(p1+p2,kq)/s;
  double yb = 2.0 * dotp(p1+p2,kqb)/s;
  double yhelp0 = 1.0/(1.0 - y - yb);
  double yhelp1 = 1.0/(1.0/2.0*s - dotp(p1,p3) - dotp(p1,kq));
  double yhelp2 = 1.0/(s - 1.0/2.0*yb*s - 1.0/2.0*y*s - dotp(p1,p3));
  double yhelp3 = 1.0/( - 1.0/2.0*s + 1.0/2.0*yb*s 
			+ dotp(p1,p3) + dotp(p1,kq));
  double yhelp4 = 1.0/(1.0/2.0*y*s - dotp(p1,kq));

  IntType MapleGenVar1,MapleGenVar2,MapleGenVar3,MapleGenVar4,MapleGenVar5,
    MapleGenVar6,MapleGenVar7;

  
#include "remainder_gg.dec"  
#include "remainder_gg.cpp"

}

