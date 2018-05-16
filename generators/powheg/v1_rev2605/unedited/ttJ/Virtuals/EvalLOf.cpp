// $Modified: Thu Oct  5 11:04:44 2006 by puwer $
#include "FourMomentum.h" 
#include "Evalf.h"
#include <cmath>

void EvalLOfgg(double f[7][133], 
	     const FourMomentum & p1, const FourMomentum & p2,
	     const FourMomentum & p3,
	     const FourMomentum & kq, const FourMomentum & kqb){
  const double m = sqrt(dotp(kq,kq));
#include "LO_gg.dec"
#include "LO_gg.cpp"  
}

void EvalLOfqq(double f[5][SMEMAXQQ], 
	     const FourMomentum & p1, const FourMomentum & p2,
	     const FourMomentum & p3,
	     const FourMomentum & kq, const FourMomentum & kqb){
  const double m = sqrt(dotp(kq,kq));
  const double s = 2.0 * dotp(p1,p2);
#include "LO_qq.dec"
#include "LO_qq.cpp"  
}
