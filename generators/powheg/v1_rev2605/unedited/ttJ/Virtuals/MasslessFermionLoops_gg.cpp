// $Modified: Thu Sep 21 12:06:28 2006 by puwer $
#include "Reduction.h"
#include "Evalf.h"
#include <complex>
#include <iostream>

using namespace std;

void MasslessFermionLoops_gg(IntType f_nfl[12][133], 
			  const FourMomentum & p1, const FourMomentum & p2, 
			  const FourMomentum & p3,
			  const FourMomentum & kq, const FourMomentum & kqb,
			  IntType* Bptr[], IntType* Cptr[],
			  IntType* Dptr[], double rterms 
			  ) {

  double nfl=1.0;

  double m = kq.getmass();

  IntType MapleGenVar1,MapleGenVar2,MapleGenVar3,MapleGenVar4,MapleGenVar5;

#include "NLOnfl_gg.dec"  

#include "NLOnfl_gg.cpp"

}

