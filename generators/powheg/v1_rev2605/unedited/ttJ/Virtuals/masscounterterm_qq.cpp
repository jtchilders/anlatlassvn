// $Modified: Thu Oct  5 11:06:59 2006 by puwer $
#include "Reduction.h"
#include "Evalf.h"
#include "ScalarInt.h"
#include <complex>
#include <iostream>

using namespace std;

double den(double);


void masscounterterm_qq(IntType f[5][SMEMAXQQ], 
		    const FourMomentum & p1, const FourMomentum & p2, 
		    const FourMomentum & p3,
		    const FourMomentum & kq, const FourMomentum & kqb,
		    double rterms
		    ) {

  const double m = kq.getmass();
  const double s = 2.0 * dotp(p1,p2);
  const double y  = 2.0 * dotp(p1+p2,kq)/s;
  const double yb = 2.0 * dotp(p1+p2,kqb)/s;
  const double yhelp0 = 1.0/(1.0 - y - yb);
  const double yhelp1 = 1.0/(1.0/2.0*s - dotp(p1,p3) - dotp(p1,kq));
  const double yhelp2 = 1.0/(s - 1.0/2.0*yb*s - 1.0/2.0*y*s - dotp(p1,p3));
  const double yhelp3 = 1.0/( - 1.0/2.0*s + 1.0/2.0*yb*s 
			+ dotp(p1,p3) + dotp(p1,kq));
  const double yhelp4 = 1.0/(1.0/2.0*y*s - dotp(p1,kq));


#include "massct_qq.dec"  
#include "massct_qq.cpp"

}

