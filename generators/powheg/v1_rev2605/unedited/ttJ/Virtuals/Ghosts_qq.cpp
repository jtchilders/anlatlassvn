// $Modified: Thu Oct  5 11:02:44 2006 by puwer $
#include "Reduction.h"
#include "Evalf.h"
#include "ScalarInt.h"
#include <complex>
#include <iostream>

using namespace std;


double den(double);

void Ghosts_qq(IntType f[5][SMEMAXQQ], 
	    const FourMomentum & p1, const FourMomentum & p2, 
	    const FourMomentum & p3,
	    const FourMomentum & kq, const FourMomentum & kqb,
	    IntType* Bptr[], IntType* Cptr[],
	    IntType* Dptr[],
	    double rterms
	    ) {

  double m = kq.getmass();
  double s = 2.0 * dotp(p1,p2);


#include "ghosts_qq.dec"  
#include "ghosts_qq.cpp"

}

