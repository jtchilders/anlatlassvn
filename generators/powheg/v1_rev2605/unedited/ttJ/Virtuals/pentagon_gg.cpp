// $Modified: Tue Jan 15 14:09:25 2008 by uwer $
#include "GGReduction.h"
#include "Reduction.h"
#include "Evalf.h"
#include <complex>
#include <iostream>

using namespace std;


static int what = 0;

void selectPentagon(int i){
  what = i;
}

#define HOHO(x_) cout << #x_ << " " << x_ << endl

void pentagon_gg(IntType f_nfl[12][133], 
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

  Topology *ETOPO11, *ETOPO12, *ETOPO13, *ETOPO14, *ETOPO15, *ETOPO16; 
  Topology *ETOPO21, *ETOPO22, *ETOPO23, *ETOPO24, *ETOPO25, *ETOPO26; 
  Topology *ETOPO31, *ETOPO32, *ETOPO33, *ETOPO34, *ETOPO35, *ETOPO36; 
  Topology *ETOPO41, *ETOPO42, *ETOPO43, *ETOPO44, *ETOPO45, *ETOPO46;

    // tij = (pi-pj)^2
  // sij = (pi+pj)^2
  const double mtq = kq.getmass2();
  const double t13 = -2.0 * dotp(p1,p3);
  const double t23 = -2.0 * dotp(p2,p3);
  const double s12 =  2.0 * dotp(p1,p2);
  const double sqqb = 2.0 * mtq + 2.0 * dotp(kq,kqb);

  const double t1q  = mtq - 2.0 * dotp(p1,kq);
  const double t1qb = mtq - 2.0 * dotp(p1,kqb);
  const double t2q  = mtq - 2.0 * dotp(p2,kq);
  const double t2qb = mtq - 2.0 * dotp(p2,kqb);
  const double s3q  = mtq + 2.0 * dotp(p3,kq);
  const double s3qb = mtq + 2.0 * dotp(p3,kqb);
  
  pentagon(ETOPO11, - 2*mtq, - 2*mtq,t13 - 2*mtq,t2q - mtq,s12 - 2*mtq, 
	   - 2*mtq,s3qb - mtq,sqqb - 2*mtq,0,0,mtq,mtq,mtq,mtq,0);

  pentagon(ETOPO12, - 2*mtq, - 2*mtq,t13 - 2*mtq,t2qb - mtq,s12 - 2*mtq, 
	   - 2*mtq,s3q - mtq,sqqb - 2*mtq,0,0,mtq,mtq,mtq,mtq,0);
  pentagon(ETOPO13, - 2*mtq, - 2*mtq,t23 - 2*mtq,t1q - mtq,s12 - 2*mtq,
	   sqqb - 2*mtq,0, - 2*mtq,s3qb - mtq,0,mtq,mtq,mtq,mtq,0);
  pentagon(ETOPO14, - 2*mtq, - 2*mtq,t23 - 2*mtq,t1qb - mtq,s12 - 2*mtq,
	   sqqb - 2*mtq,0, - 2*mtq,s3q - mtq,0,mtq,mtq,mtq,mtq,0);
  pentagon(ETOPO15, - 2*mtq, - 2*mtq,t23 - 2*mtq,t1q - mtq,t13 - 2*mtq,
	   sqqb - 2*mtq,0, - 2*mtq,t2qb - mtq,0,mtq,mtq,mtq,mtq,0);
  pentagon(ETOPO16, - 2*mtq, - 2*mtq,t23 - 2*mtq,t1qb - mtq,t13 - 2*mtq,
	   sqqb - 2*mtq,0, - 2*mtq,t2q - mtq,0,mtq,mtq,mtq,mtq,0);
  pentagon(ETOPO21, - 2*mtq, - 2*mtq,t1q - mtq,t2qb - mtq,s12 - 2*mtq,0,
	   s3q - mtq,s3qb - mtq,0,0,mtq,mtq,mtq,0,0);
  pentagon(ETOPO22, - 2*mtq, - 2*mtq,t1qb - mtq,t2q - mtq,s12 - 2*mtq,0,
	   s3qb - mtq,s3q - mtq,0,0,mtq,mtq,mtq,0,0);
  pentagon(ETOPO23, - 2*mtq, - 2*mtq,t1q - mtq,s3qb - mtq,t13 - 2*mtq,0,
	   t2q - mtq,t2qb - mtq,0,0,mtq,mtq,mtq,0,0);
  pentagon(ETOPO24, - 2*mtq, - 2*mtq,t1qb - mtq,s3q - mtq,t13 - 2*mtq,0,
	   t2qb - mtq,t2q - mtq,0,0,mtq,mtq,mtq,0,0);
  pentagon(ETOPO25,0,t1qb - mtq,t2q - mtq,0,t23 - 2*mtq, - 2*mtq,t1q - mtq, 
	   - 2*mtq,0,s3qb - mtq,0,mtq,mtq,mtq,0);
  pentagon(ETOPO26,0,t1qb - mtq,s3q - mtq,0,t23 - 2*mtq, - 2*mtq,t1q - mtq, 
	   - 2*mtq,0,t2qb - mtq,0,mtq,mtq,mtq,0);
  pentagon(ETOPO31,t1qb - mtq,t2q - mtq,0,0, - 2*mtq,0,s3q - mtq,s3qb - mtq,
	   0,s12,0,mtq,mtq,0,0);
  pentagon(ETOPO32,t1qb - mtq,s3q - mtq,0,0, - 2*mtq,0,t2q - mtq,t2qb - mtq,
	   0,t13,0,mtq,mtq,0,0);
  pentagon(ETOPO33, - 2*mtq,0,t1qb - mtq,t2q - mtq,t1q - mtq,0,s3qb - mtq,
	   t23,0,0,mtq,mtq,0,0,0);
  pentagon(ETOPO34, - 2*mtq,0,t1qb - mtq,s3q - mtq,t1q - mtq,0,t2qb - mtq,
	   t23,0,0,mtq,mtq,0,0,0);
  pentagon(ETOPO35,t1q - mtq,t2qb - mtq,0,0, - 2*mtq,0,s3qb - mtq,s3q - mtq,
	   0,s12,0,mtq,mtq,0,0);
  pentagon(ETOPO36,t1q - mtq,s3qb - mtq,0,0, - 2*mtq,0,t2qb - mtq,t2q - mtq,
	   0,t13,0,mtq,mtq,0,0);
  pentagon(ETOPO41,t1q - mtq,0,0,t23,0,s3qb - mtq,0,s12,sqqb,0,0,mtq,0,0,0);
  pentagon(ETOPO42,t1q - mtq,0,0,t23,0,t2qb - mtq,0,t13,sqqb,0,0,mtq,0,0,0);
  pentagon(ETOPO43,t1qb - mtq,0,0,t23,0,s3q - mtq,0,s12,sqqb,0,0,mtq,0,0,0);
  pentagon(ETOPO44,t2q - mtq,0,0,t13,s3qb - mtq,0,0,s12,0,sqqb,0,mtq,0,0,0);
  pentagon(ETOPO45,t1qb - mtq,0,0,t23,0,t2q - mtq,0,t13,sqqb,0,0,mtq,0,0,0);
  pentagon(ETOPO46,t2qb - mtq,0,0,t13,s3q - mtq,0,0,s12,0,sqqb,0,mtq,0,0,0);

#ifdef ORIVERSION
  FourMomentum O("0",0.0);
  O.setFourMomentum(0.0,0.0,0.0,0.0);

  ETOPO11 = new Topology(pentagon(O,-p1,p2, - p1 + p3,p2-kq,m,m,m,m,0));
  ETOPO12 = new Topology(pentagon(O,-p1,p2, - p1 + p3,p2 - kqb,m,m,m,m,0));
  ETOPO13 = new Topology(pentagon(O,-p1,p2,p2 - p3, - p1 + kq,m,m,m,m,0));
  ETOPO14 = new Topology(pentagon(O,-p1,p2,p2 - p3, - p1 + kqb,m,m,m,m,0));
  ETOPO15 = new Topology(pentagon(O,-p1,-p3,p2 - p3, - p1 + kq,m,m,m,m,0));  
  ETOPO16 = new Topology(pentagon(O,-p1,-p3,p2 - p3, - p1 + kqb,m,m,m,m,0));

  ETOPO21 = new Topology(pentagon(O,-p1,p2, - p1 + kq,p2 - kqb,m,m,m,0,0));
  ETOPO22 = new Topology(pentagon(O,-p1,p2, - p1 + kqb,p2 - kq,m,m,m,0,0));
  ETOPO23 = new Topology(pentagon(O,-p1,-p3, - p1 + kq, - p3 - kqb,m,m,m,0,0));
  ETOPO24 = new Topology(pentagon(O,-p1,-p3, - p1 + kqb, - p3 - kq,m,m,m,0,0));
  ETOPO25 = new Topology(pentagon(O,-kq, - p1 + kqb,p2 - kq,-p1,0,m,m,m,0));
  ETOPO26 = new Topology(pentagon(O,-kq, - p1 + kqb, - p3 - kq,-p1,0,m,m,m,0));
  
  ETOPO31 = new Topology(pentagon(O, - p1 + kqb,p2 - kq,-p1,p2,0,m,m,0,0));
  ETOPO32 = new Topology(pentagon(O, - p1 + kqb, - p3 - kq,-p1,-p3,0,m,m,0,0));
  ETOPO33 = new Topology(pentagon(O,-p1,-kq, - p1 + kqb,p2 - kq,m,m,0,0,0));
  ETOPO34 = new Topology(pentagon(O,-p1,-kq, - p1 + kqb, - p3 - kq,m,m,0,0,0));
  ETOPO35 = new Topology(pentagon(O, - p1 + kq,p2 - kqb,-p1,p2,0,m,m,0,0));
  ETOPO36 = new Topology(pentagon(O, - p1 + kq, - p3 - kqb,-p1,-p3,0,m,m,0,0));
  
  ETOPO41 = new Topology(pentagon(O, - p1 + kq,-p1,p2,p2 - p3,0,m,0,0,0));
  ETOPO42 = new Topology(pentagon(O, - p1 + kq,-p1,-p3,p2 - p3,0,m,0,0,0));
  ETOPO43 = new Topology(pentagon(O, - p1 + kqb,-p1,p2,p2 - p3,0,m,0,0,0));
  ETOPO44 = new Topology(pentagon(O,p2 - kq,-p1,p2, - p1 + p3,0,m,0,0,0));
  ETOPO45 = new Topology(pentagon(O, - p1 + kqb,-p1,-p3,p2 - p3,0,m,0,0,0));
  ETOPO46 = new Topology(pentagon(O,p2 - kqb,-p1,p2, - p1 + p3,0,m,0,0,0));
#endif

  Topology 
    *ETOPO11x1=ETOPO11->subtopo(1),
    *ETOPO12x1=ETOPO12->subtopo(1),
    *ETOPO13x1=ETOPO13->subtopo(1),
    *ETOPO14x1=ETOPO14->subtopo(1),
    *ETOPO15x1=ETOPO15->subtopo(1),
    *ETOPO16x1=ETOPO16->subtopo(1),
    
    *ETOPO21x1=ETOPO21->subtopo(1),
    *ETOPO22x1=ETOPO22->subtopo(1),
    *ETOPO23x1=ETOPO23->subtopo(1),
    *ETOPO24x1=ETOPO24->subtopo(1),
    *ETOPO25x1=ETOPO25->subtopo(1),
    *ETOPO26x1=ETOPO26->subtopo(1),
    
    *ETOPO31x1=ETOPO31->subtopo(1),
    *ETOPO32x1=ETOPO32->subtopo(1),
    *ETOPO33x1=ETOPO33->subtopo(1),
    *ETOPO34x1=ETOPO34->subtopo(1),
    *ETOPO35x1=ETOPO35->subtopo(1),
    *ETOPO36x1=ETOPO36->subtopo(1),
  
    *ETOPO41x1=ETOPO41->subtopo(1),
    *ETOPO42x1=ETOPO42->subtopo(1),
    *ETOPO43x1=ETOPO43->subtopo(1),
    *ETOPO44x1=ETOPO44->subtopo(1),
    *ETOPO45x1=ETOPO45->subtopo(1),
    *ETOPO46x1=ETOPO46->subtopo(1);


  switch (what){
  case ALLPENTAGONS:
    // Keep all the topologies
    break;

  case PENTA1:
    {
      ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
      ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;
      
      ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
      ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;
      
      ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
      ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
      cout << "Only pentagons with one internal gluon line are kept\n";
    }
    break;
  case PENTA11:
    ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;
    
    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;
    
    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;
    
    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO11->getDetS() << endl;
    break;
  case PENTA12:
    ETOPO11 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;
    
    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;
    
    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;
    
    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO12->getDetS() << endl;
    break;
  case PENTA13:
    ETOPO11 = ETOPO12 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;
    
    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;
    
    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;
    
    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO13->getDetS() << endl;
    break;
  case PENTA14:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO15x1 = ETOPO16x1 = 0;
    
    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;
    
    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;
    
    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO14->getDetS() << endl;
    break;
  case PENTA15:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1  = ETOPO16x1 = 0;
    
    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;

    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;
 
    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO15->getDetS() << endl;
     break;
  case PENTA16:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = 0;
    
    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;
    
    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;
    
    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO16->getDetS() << endl;
    break;
  case PENTA2: 
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;

    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;
    
    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    break;
  case PENTA21:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;
 
    ETOPO22 = ETOPO23  = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;

    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;
    
    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO21->getDetS() << endl;
    break;  
  case PENTA22:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;
 
    ETOPO21 = ETOPO23  = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;

    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;
    
    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO22->getDetS() << endl;
    break;
  case PENTA23:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;
 
    ETOPO21 = ETOPO22  = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;

    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;
    
    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO23->getDetS() << endl;
    break;
  case PENTA24:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;
 
    ETOPO21 = ETOPO22  = ETOPO23 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO25x1 = ETOPO26x1 = 0;

    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;
    
    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO24->getDetS() << endl;
    break;
  case PENTA25:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;
 
    ETOPO21 = ETOPO22  = ETOPO23 = ETOPO24 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO26x1 = 0;

    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;
    
    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO25->getDetS() << endl;
    break;
  case PENTA26:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;
 
    ETOPO21 = ETOPO22  = ETOPO23 = ETOPO24 = ETOPO25 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = 0;

    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;
    
    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO26->getDetS() << endl;
    break;
  case PENTA3:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;

    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;
 
    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    break;
  case PENTA31:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;

    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;

    ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0; 

    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO31->getDetS() << endl;
    break;  
case PENTA32:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;

    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;

    ETOPO31 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0; 

    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO32->getDetS() << endl;
    break;
case PENTA33:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;

    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;

    ETOPO31 = ETOPO32 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0; 

    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO33->getDetS() << endl;
    break;
case PENTA34:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;

    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;

    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO35x1 = ETOPO36x1 = 0; 

    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO34->getDetS() << endl;
    break;
case PENTA35:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;

    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;

    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO36x1 = 0; 

    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO35->getDetS() << endl;
    break;
  case PENTA36:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;

    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;

    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = 0; 

    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO36->getDetS() << endl;
    break;
  
  case PENTA4:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;

    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;
 
    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;
    break;
  case PENTA41:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;

    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;
 
    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;
   
    ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO41->getDetS() << endl;
    break;  
  case PENTA42:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;

    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;
 
    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;

    ETOPO41 = ETOPO43 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO42->getDetS() << endl;
    break;  

  case PENTA43:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;

    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;
 
    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;

    ETOPO41 = ETOPO42 = ETOPO44 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO44x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO43->getDetS() << endl;
    break;  

  case PENTA44:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;

    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;
 
    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;

    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO45 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO45x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO44->getDetS() << endl;
    break;  

  case PENTA45:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;

    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;
 
    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;

    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO46 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO46x1 = 0;
    cout << "detS() = " << ETOPO45->getDetS() << endl;
    break;

  case PENTA46:
    ETOPO11 = ETOPO12 = ETOPO13 = ETOPO14 = ETOPO15 = ETOPO16 = 0;
    ETOPO11x1 = ETOPO12x1 = ETOPO13x1 = ETOPO14x1 = ETOPO15x1 = ETOPO16x1 = 0;

    ETOPO21 = ETOPO22 = ETOPO23 = ETOPO24 = ETOPO25 = ETOPO26 = 0;
    ETOPO21x1 = ETOPO22x1 = ETOPO23x1 = ETOPO24x1 = ETOPO25x1 = ETOPO26x1 = 0;
 
    ETOPO31 = ETOPO32 = ETOPO33 = ETOPO34 = ETOPO35 = ETOPO36 = 0;
    ETOPO31x1 = ETOPO32x1 = ETOPO33x1 = ETOPO34x1 = ETOPO35x1 = ETOPO36x1 = 0;

    ETOPO41 = ETOPO42 = ETOPO43 = ETOPO44 = ETOPO45 = 0;
    ETOPO41x1 = ETOPO42x1 = ETOPO43x1 = ETOPO44x1 = ETOPO45x1 = 0;
    cout << "detS() = " << ETOPO46->getDetS() << endl;
    break;
  }


#include "penta_gg.dec"  
#include "penta_gg.cpp"

  Topology::cacheclear();
  /*  

      delete ETOPO11;
  delete ETOPO12; 
  delete ETOPO13;
  delete ETOPO14;
  delete ETOPO15;
  delete ETOPO16;
  delete ETOPO21; 
  delete ETOPO22; 
  delete ETOPO23; 
  delete ETOPO24; 
  delete ETOPO25; 
  delete ETOPO26; 
  delete ETOPO31; 
  delete ETOPO32; 
  delete ETOPO33; 
  delete ETOPO34; 
  delete ETOPO35; 
  delete ETOPO36; 
  delete ETOPO41; 
  delete ETOPO42; 
  delete ETOPO43; 
  delete ETOPO44; 
  delete ETOPO45;
  delete ETOPO46; 
  */
  /*  delete ETOPO11x1, ETOPO12x1,  ETOPO13x1, ETOPO14x1, 
    ETOPO21x1, ETOPO22x1, ETOPO23x1, ETOPO24x1, ETOPO25x1, ETOPO26x1, 
    ETOPO31x1, ETOPO32x1, ETOPO33x1, ETOPO34x1, ETOPO35x1, ETOPO36x1, 
    ETOPO41x1, ETOPO42x1, ETOPO43x1, ETOPO44x1, ETOPO45x1, ETOPO46x1;*/

}

