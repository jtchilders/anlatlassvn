// $Modified: Thu Oct  5 11:06:16 2006 by puwer $
// g++ -c -I~/lib/c++/Tools EvalSMExx.cpp
#include "FourMomentum.h"
#include "Evalf.h"
#include "EvalSMExx.h"


using namespace std;

void EvalSMEgg(complex<double>S[2][2][2][133],
	     const FourMomentum & p1, const FourMomentum & p2, 
	     const FourMomentum & p3, const FourMomentum & q1,
	     const FourMomentum & q2, const FourMomentum & r1, 
	     const FourMomentum & r2){
  double m = sqrt(2.0*dotp(q1,q2));
#include "SME_gg.dec"
#include "SME_gg.cpp"      
}

void EvalSMEqq(complex<double>S[2][2][SMEMAXQQ],
	     const FourMomentum & p1, const FourMomentum & p2, 
	     const FourMomentum & p3, const FourMomentum & q1,
	     const FourMomentum & q2, const FourMomentum & r1, 
	     const FourMomentum & r2){
  const double sqrt2 = sqrt(2.);
  double m = sqrt(2.0*dotp(q1,q2));
#include "SME_qq.dec"
#include "SME_qq.cpp"      
}
