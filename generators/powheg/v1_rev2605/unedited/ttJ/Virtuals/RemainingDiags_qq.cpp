// $Modified: Thu Oct  5 11:02:57 2006 by puwer $
#include "Reduction.h"
#include "Evalf.h"
#include "ScalarInt.h"
#include <complex>
#include <iostream>

using namespace std;



void RemainingDiags_qq(IntType f[5][SMEMAXQQ], 
		    const FourMomentum & p1, const FourMomentum & p2, 
		    const FourMomentum & p3,
		    const FourMomentum & kq, const FourMomentum & kqb,
		    IntType* Bptr[], IntType* Cptr[],
		    IntType* Dptr[],
		    double rterms
		    ) {

  double s = 2.0 * dotp(p1,p2);
  double m = kq.getmass();

  IntType MapleGenVar1,MapleGenVar2,MapleGenVar3,MapleGenVar4,MapleGenVar5,
    MapleGenVar6,MapleGenVar7,MapleGenVar8,MapleGenVar9;

#include "remainder_qq.dec"  
#include "remainder_qq.cpp"

}

