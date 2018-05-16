// $Modified: Thu Oct  5 11:02:51 2006 by puwer $
#include "Reduction.h"
#include "Evalf.h"
#include <complex>
#include <iostream>

using namespace std;

void MasslessFermionLoops_qq(IntType f[5][SMEMAXQQ], 
			  const FourMomentum & p1, const FourMomentum & p2, 
			  const FourMomentum & p3,
			  const FourMomentum & kq, const FourMomentum & kqb,
			  IntType* Bptr[], IntType* Cptr[],
			  IntType* Dptr[], double rterms 
			  ) {

  double nfl=1.0;

  double m = kq.getmass();
  double s = 2.0 * dotp(p1,p2); 
  IntType MapleGenVar1,MapleGenVar2,MapleGenVar3,MapleGenVar4,MapleGenVar5;

#include "NLOnfl_qq.dec"  

#include "NLOnfl_qq.cpp"

}

