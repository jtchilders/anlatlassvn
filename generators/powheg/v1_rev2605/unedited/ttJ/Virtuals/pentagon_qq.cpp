// $Modified: Tue Jan 15 14:09:38 2008 by uwer $
#include "GGReduction.h"
#include "Reduction.h"
#include "Evalf.h"
#include <complex>
#include <iostream>

using namespace std;



#define HOHO(x_) cout << #x_ << " " << x_ << endl

void pentagon_qq(IntType f[5][SMEMAXQQ], 
	  const FourMomentum & p1, const FourMomentum & p2, 
	  const FourMomentum & p3,
	  const FourMomentum & kq, const FourMomentum & kqb,
	  double rterms) {

  double m = kq.getmass();
  

  double s = 2.0 * dotp(p1,p2);
  double y = 2.0 * dotp(p1+p2,kq)/s;
  double yb = 2.0 * dotp(p1+p2,kqb)/s;
  double MapleGenVar1,MapleGenVar2,MapleGenVar3,MapleGenVar4,MapleGenVar5,
    MapleGenVar6,MapleGenVar7;

  Topology *ETOPO1, *ETOPO2, *ETOPO3, *ETOPO4, *ETOPO5, *ETOPO6,
    *ETOPO7, *ETOPO8; 

    // tij = (pi-pj)^2
  // sij = (pi+pj)^2

  const double mtq = kq.getmass2();
  const double xt13 = -2.0 * dotp(p1,p3);
  const double xt23 = -2.0 * dotp(p2,p3);
  const double xs12 =  2.0 * dotp(p1,p2);
  const double xsqqb = 2.0 * mtq + 2.0 * dotp(kq,kqb);

  const double xt1q  = mtq - 2.0 * dotp(p1,kq);
  const double xt1qb = mtq - 2.0 * dotp(p1,kqb);
  const double xt2q  = mtq - 2.0 * dotp(p2,kq);
  const double xt2qb = mtq - 2.0 * dotp(p2,kqb);
  const double xs3q  = mtq + 2.0 * dotp(p3,kq);
  const double xs3qb = mtq + 2.0 * dotp(p3,kqb);
  
  pentagon(ETOPO1,xt1q - mtq,xt2qb - mtq,0,0, - 2*mtq,0,xs3qb - mtq,
	   xs3q - mtq,0,xs12,0,mtq,mtq,0,0);

  pentagon(ETOPO2,xt1qb - mtq,xt2q - mtq,0,0, - 2*mtq,0,xs3q - mtq,
	   xs3qb - mtq,0,xs12,0,mtq,mtq,0,0);
  
  pentagon(ETOPO3,xt1q - mtq,0,xt23,0,0,0,xt2qb - mtq,xsqqb,
	   xt13,0,0,mtq,0,0,0);

  pentagon(ETOPO4,xt1qb - mtq,0,xt23,0,0,0,xt2q - mtq,xsqqb,xt13,0,0,mtq,0,0,0);

  pentagon(ETOPO5,xt1q - mtq,0,0,xt23,0,xs3qb - mtq,0,xs12,xsqqb,
	   0,0,mtq,0,0,0);

  pentagon(ETOPO6,xt1qb - mtq,0,0,xt23,0,xs3q - mtq,0,xs12,xsqqb,
	   0,0,mtq,0,0,0);

  pentagon(ETOPO7,xt2q - mtq,0,0,xt13,xs3qb - mtq,0,0,xs12,0,xsqqb,
	   0,mtq,0,0,0);

  pentagon(ETOPO8,xt2qb - mtq,0,0,xt13,xs3q - mtq,0,0,xs12,0,xsqqb,
	   0,mtq,0,0,0);


#include "penta_qq.dec"  
#include "penta_qq.cpp"

  Topology::cacheclear();

}

