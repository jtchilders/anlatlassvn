// $Modified: Mon Aug  6 17:38:48 2007 by uwer $
#include "StandardModelParameters.h"
#include "EvalSMExx.h"
#include "FourMomentum.h"
#include "Subtractions.h"
#include "Color.h"
#include "Reduction.h"
#include "GGReduction.h"
#include "ScalarInt.h"
#include "Errors.h"
#include "Evalf.h"
#include <complex>
#include <cmath>
#include <cstring>
#include <iostream>
#include "EPSHisto.h"

using namespace std;


static StandardModelParameters& parms = StandardModelParameters::instance();



class KahanSum {
private:
  typedef complex<double> itype;
  itype ct,cy,cc,csum;
public:
  KahanSum(){
    reset();
  }
  void reset(){
    csum = itype(0.0,0.0);
    cc = itype(0.0,0.0); 
  }
  void add(complex<double> x){
    cy = static_cast<itype>(x) - cc;
    ct = csum + cy;
    cc = (ct - csum); 
    cc -= cy;
    csum = ct;
  }
  complex<double> getSum(){
    return( static_cast< complex<double> >(csum) );
  }
};

void setNLOtozero(IntType f[12][133]) {
  for(int color=1; color < 12; color++){
    for(int i=1;i<133;i++){
      f[color][i] = 0.0;
    }
  }
}

double gggtt2(const FourMomentum & p1, const FourMomentum & p2, 
	      const FourMomentum & p3,
	      const FourMomentum & kq, const FourMomentum & kqb){

  FourMomentum q1, q2, r1, r2; 

  complex<double>S[2][2][2][133];

  double f[7][133];  
  complex<double> ampLO[7], tmp;

  double amp2=0;

  static double CLO[7][7];
  static bool initialized=false;


  if ( !initialized ) {
    LOcolorgg(CLO);
    initialized = true;
  }

  evalqr(kq,kqb,q1,q2,r1,r2);
  EvalLOfgg(f, p1,p2,p3,kq,kqb);



  static KahanSum summer;
  /*
   *  Sum over the polarization and square it:
   */

  EvalSMEgg(S,p1, p2, p3,q1,q2,r1,r2);

  for (int tspin=0; tspin < 2; tspin++) {

    if (tspin == 1) 
      EvalSMEgg(S,p1, p2, p3,q1,q2,r2,r1);

    if (tspin == 2) 
      EvalSMEgg(S,p1, p2, p3,q2,q1,r1,r2);

    if (tspin == 3) 
      EvalSMEgg(S,p1, p2, p3,q2,q1,r2,r1);

    for (int pol1=0; pol1 < 2; pol1++){
      for (int pol2=0; pol2 < 2; pol2++){
	for (int pol3=0; pol3 < 2; pol3++){
	  for(int color=1;color < 7;color++){
	    summer.reset();
	    for(int i=1; i < 133; i++){
	      // ampLO[color] +=  S[pol1][pol2][pol3][i] * f[color][i];
	      summer.add( S[pol1][pol2][pol3][i] * f[color][i] );
	    }
	    ampLO[color] = summer.getSum();
	  }
	  summer.reset();
	  for(int c1=1; c1 < 7; c1++) {
	    for(int c2=1; c2 < 7; c2++) {
	      //tmp = tmp + conj(ampLO[c1]) * CLO[c1][c2] * ampLO[c2];
	      summer.add(conj(ampLO[c1]) * CLO[c1][c2] * ampLO[c2]);
	    }
	  }
	  //	  cout << tmp << endl;
	  amp2 = amp2 + summer.getSum().real();
	}
      }
    }
  }
  return 2.0*amp2;
}

void iopgg(complex<double> Mij[7][7],
	   const FourMomentum & p1, const FourMomentum & p2, 
	   const FourMomentum & p3,
	   const FourMomentum & kq, const FourMomentum & kqb ) {

  double m = kq.getmass();

  static double clo[7][7];
  static double cg1g2[7][7];
  static double cqqb[7][7];
  static double cg1g3[7][7]; 
  static double cg2g3[7][7]; 
  static double cg1q[7][7];
  static double cg2q[7][7]; 
  static double cg3q[7][7]; 
  static double cg1qb[7][7];
  static double cg2qb[7][7]; 
  static double cg3qb[7][7];

  static bool initialized = false;

  double mu2 = parms.getRenormalizationScaleSquared();

  if (!initialized) {
    ggcorrelations(clo, cg1g2, cqqb,cg1g3, cg2g3, cg1q,
		   cg2q, cg3q, cg1qb,cg2qb, cg3qb);
    initialized = true;
  }


  const double kappa = (2.0/3.0);
  vector<double> mf;

  for(int c1=1; c1 < 7; c1++) {
    for(int c2=c1; c2 < 7; c2++) {

      using namespace Subtractions;
  
      Mij[c1][c2] = 
	-1.0 *  (
	   // Final state I_m
	   + cqqb[c1][c2] 
	   * ( Vq( 2.0*dotp(kq,kqb),m,m) - (pi2/3.0) 
	       + ( Gamma_q(m) + gammaq*log(mu2/2.0/dotp(kq,kqb)) 
		   + gammaq + Kq )/CF 
	       )
	   
	   + cqqb[c1][c2] 
	   * ( Vq( 2.0*dotp(kq,kqb),m,m) - (pi2/3.0) 
	       + ( Gamma_q(m) + gammaq*log(mu2/2.0/dotp(kq,kqb)) 
		   + gammaq + Kq )/CF 
	       )
	   
	   + cg3q[c1][c2] 
	   * ( Vq( 2.0*dotp(kq,p3),m,0.) - (pi2/3.0) 
	       + ( Gamma_q(m) + gammaq*log(mu2/2.0/dotp(kq,p3)) 
		   + gammaq + Kq)/CF 
	       ) 
	   
	   + cg3qb[c1][c2] 
	   * ( Vq( 2.0*dotp(kqb,p3),m,0.) - (pi2/3.0) 
	       + ( Gamma_q(m) + gammaq*log(mu2/2.0/dotp(kqb,p3)) 
		   + gammaq + Kq)/CF 
	       )
	   
	   + cg3q[c1][c2] 
	   * ( Vg( 2.0*dotp(p3,kq), 0., m, mf, kappa) - (pi2/3.0) 
	       + ( Gamma_g(mf) + gammag*log(mu2/2.0/dotp(p3,kq)) 
		   + gammag + Kg )/CA 
	       )
	   
	   + cg3qb[c1][c2] 
	   * ( Vg( 2.0*dotp(p3,kqb), 0., m, mf, kappa) - (pi2/3.0) 
	       + ( Gamma_g(mf) + gammag*log(mu2/2.0/dotp(p3,kqb)) 
		   + gammag + Kg )/CA 
	       )
	   
	   // initial -- final I_a 
	   /*
	    * Final state emitter, initial state spectator
	    */
	   + cg1q[c1][c2] 
	   * ( Vq( 2.0*dotp(kq,p1),m,0.) - (pi2/3.0) 
	       + ( Gamma_q(m) + gammaq*log(mu2/2.0/dotp(kq,p1)) 
		   + gammaq + Kq )/CF 
	       )
	   
	   + cg1qb[c1][c2] 
	   * ( Vq( 2.0*dotp(kqb,p1),m,0.) - (pi2/3.0) 
	       + ( Gamma_q(m) + gammaq*log(mu2/2.0/dotp(kqb,p1)) 
		   + gammaq + Kq )/CF 
	       )
	   
	   + cg1g3[c1][c2] 
	   * ( Vg( 2.0*dotp(p3,p1), 0., 0., mf, kappa) - (pi2/3.0) 
	       + ( Gamma_g(mf) + gammag*log(mu2/2.0/dotp(p3,p1)) 
		   + gammag + Kg )/CA 
	       )

	   /*
	    * Final state spectator, initial state emitter
	    */
	   + cg1q[c1][c2] 
	   * ( Vg( 2.0*dotp(p1,kq), 0., m, mf, 2.0/3.0) - (pi2/3.0) 
	       + ( Gamma_g() + gammag*log(mu2/2.0/dotp(p1,kq)) 
		   + gammag + Kg )/CA 
	       )
	   
	   + cg1qb[c1][c2] 
	   * ( Vg( 2.0*dotp(p1,kqb), 0., m, mf, 2.0/3.0) - (pi2/3.0) 
	       + ( Gamma_g() + gammag*log(mu2/2.0/dotp(p1,kqb)) 
		   + gammag + Kg )/CA 
	       )
	   + cg1g3[c1][c2] 
	   * ( Vg( 2.0*dotp(p1,p3), 0., 0., mf, 2.0/3.0) - (pi2/3.0) 
	       + ( Gamma_g() + gammag*log(mu2/2.0/dotp(p1,p3)) 
		   + gammag + Kg )/CA 
	       )
	   // initial -- final I_a 
	   /*
	    * Final state emitter, initial state spectator
	    */
	   + cg2q[c1][c2] 
	   * ( Vq( 2.0*dotp(kq,p2),m,0.) - (pi2/3.0) 
	       + ( Gamma_q(m) + gammaq*log(mu2/2.0/dotp(kq,p2)) 
		   + gammaq + Kq )/CF 
	       )
	   
	   + cg2qb[c1][c2] 
	   * ( Vq( 2.0*dotp(kqb,p2),m,0.) - (pi2/3.0) 
	       + ( Gamma_q(m) + gammaq*log(mu2/2.0/dotp(kqb,p2)) 
		   + gammaq + Kq)/CF 
	       )
	   
	   + cg2g3[c1][c2] 
	   * ( Vg( 2.0*dotp(p3,p2), 0., 0., mf, kappa) - (pi2/3.0) 
	       + ( Gamma_g(mf) + gammag*log(mu2/2.0/dotp(p3,p2)) 
		   + gammag + Kg )/CA 
	       )
	   /*
	    * Final state spectator, initial state emitter
	    */
	   + cg2q[c1][c2] 
	   * ( Vg( 2.0*dotp(p2,kq), 0., m, mf, 2.0/3.0) - (pi2/3.0) 
	   + ( Gamma_g() + gammag*log(mu2/2.0/dotp(p2,kq)) 
	       + gammag + Kg )/CA 
	       )
	   
	   + cg2qb[c1][c2] 
	   * ( Vg( 2.0*dotp(p2,kqb), 0., m, mf, 2.0/3.0) - (pi2/3.0) 
	       + ( Gamma_g() + gammag*log(mu2/2.0/dotp(p2,kqb)) 
		   + gammag + Kg )/CA 
	       )
	   + cg2g3[c1][c2] 
	   * ( Vg( 2.0*dotp(p2,p3), 0., 0., mf, 2.0/3.0) - (pi2/3.0) 
	       + ( Gamma_g() + gammag*log(mu2/2.0/dotp(p2,p3)) 
		   + gammag + Kg )/CA 
	       )
	   // The two additional contributions from  
	   // [hep-ph/0201036, eq. 6.66] 
	   + cg1g2[c1][c2] 
	   * ( Vsing(2.0*dotp(p1,p2), 0., 0.) - (pi2/3.0) 
	       + ( Gamma_g() + gammag * log(mu2/2.0/dotp(p1,p2)) 
		   + gammag + Kg )/CA
	       )
	   
	   + cg1g2[c1][c2] 
	   * ( Vsing(2.0*dotp(p1,p2), 0., 0.) - (pi2/3.0) 
	       + ( Gamma_g() + gammag * log(mu2/2.0/dotp(p1,p2)) 
		   + gammag + Kg )/CA
	       )
	   );   
      Mij[c2][c1] = Mij[c1][c2];
    }
  }
}
    
double Ioperator_gg(const FourMomentum & p1, const FourMomentum & p2, 
		 const FourMomentum & p3,
		 const FourMomentum & kq, const FourMomentum & kqb){

  FourMomentum q1, q2, r1, r2; 


  complex<double>S[2][2][2][133];

  double f[7][133];

  complex<double> ampLO[7];

  complex<double> Mij[7][7];

  complex<double> tmp;

  double amp2=0;


  static bool initialized=false;


  evalqr(kq,kqb,q1,q2,r1,r2);
  EvalLOfgg(f, p1,p2,p3,kq,kqb);

  iopgg(Mij,p1,p2,p3,kq,kqb);

  /*
   *  Sum over the polarization and square it:
   */
  EvalSMEgg(S,p1, p2, p3,q1,q2,r1,r2);

  for (int tspin=0; tspin < 2; tspin++) {
    if (tspin == 1) 
      EvalSMEgg(S,p1, p2, p3,q1,q2,r2,r1);

    if (tspin == 2) 
      EvalSMEgg(S,p1, p2, p3,q2,q1,r1,r2);

    if (tspin == 3) 
      EvalSMEgg(S,p1, p2, p3,q2,q1,r2,r1);

    for (int pol1=0; pol1 < 2; pol1++){
      for (int pol2=0; pol2 < 2; pol2++){
	for (int pol3=0; pol3 < 2; pol3++){
	  for(int color=1;color < 7;color++){
	    ampLO[color] = 0.0;
	    for(int i=1; i < 133; i++){
	      ampLO[color] +=  S[pol1][pol2][pol3][i] * f[color][i];
	    }
	  }

	  tmp = complex<double>(0.,0.);
	  for(int c1=1; c1 < 7; c1++){
	    for(int c2=1; c2 < 7; c2++){
	      tmp += conj(ampLO[c1]) * Mij[c1][c2] * ampLO[c2];
	    }
	  }

	  amp2 = amp2 + tmp.real();
	}
      }
    }
  }
  return(2 * amp2);
}



static void evaluatePVcoefficients_gg(CoeffCache& PVcoeffs,
			    const FourMomentum & p1,const FourMomentum & p2,
			    const FourMomentum & p3,
			    const FourMomentum & kq, const FourMomentum & kqb,
			    IntType* Bptr[], IntType* Cptr[], IntType* Dptr[]){
  /*
   *  Calculate the PV reduction coefficients:
   */
  const double mq = kq.getmass2();  
  // tij = (pi-pj)^2
  // sij = (pi+pj)^2
  const double t13 = -2.0 * dotp(p1,p3);
  const double t23 = -2.0 * dotp(p2,p3);
  const double s12 =  2.0 * dotp(p1,p2);
  const double sqqb = 2.0 * mq + 2.0 * dotp(kq,kqb);

  const double t1q  = mq - 2.0 * dotp(p1,kq);
  const double t1qb = mq - 2.0 * dotp(p1,kqb);
  const double t2q  = mq - 2.0 * dotp(p2,kq);
  const double t2qb = mq - 2.0 * dotp(p2,kqb);
  const double s3q  = mq + 2.0 * dotp(p3,kq);
  const double s3qb = mq + 2.0 * dotp(p3,kqb);

  /*
   *  Massless fermion loops
   */
  PVcoeffs.lookup(Dptr[1],0,s12,0,t23,0,sqqb,0,0,0,0);
  PVcoeffs.lookup(Dptr[2],0,t13,0,t23,0,sqqb,0,0,0,0);
  PVcoeffs.lookup(Dptr[3],0,s12,sqqb,t13,0,0,0,0,0,0);

  PVcoeffs.lookup(Cptr[1],0,sqqb,t23,0,0,0);
  PVcoeffs.lookup(Cptr[2],0,s12,0,0,0,0);
  PVcoeffs.lookup(Cptr[3],0,t13,0,0,0,0);
  PVcoeffs.lookup(Cptr[4],0,sqqb,t13,0,0,0);
  PVcoeffs.lookup(Cptr[5],0,t23,0,0,0,0);
  PVcoeffs.lookup(Cptr[6],0,sqqb,s12,0,0,0);


  PVcoeffs.lookup(Bptr[1], t13, 0,0);
  PVcoeffs.lookup(Bptr[2], s12, 0,0);
  PVcoeffs.lookup(Bptr[3], t23, 0,0);
  PVcoeffs.lookup(Bptr[4], sqqb,0,0);

  

  /*
   *  Top quark loops
   */
  PVcoeffs.lookup(Dptr[4],0,s12,sqqb,t13,0,0,mq,mq,mq,mq);
  PVcoeffs.lookup(Dptr[5],0,t13,0,t23,0,sqqb,mq,mq,mq,mq);
  PVcoeffs.lookup(Dptr[6],0,s12,0,t23,0,sqqb,mq,mq,mq,mq);

  PVcoeffs.lookup(Cptr[7],0,sqqb,t23,mq,mq,mq);
  PVcoeffs.lookup(Cptr[8],0,s12,0,mq,mq,mq);
  PVcoeffs.lookup(Cptr[9],0,t13,0,mq,mq,mq);
  PVcoeffs.lookup(Cptr[10],0,sqqb,t13,mq,mq,mq);
  PVcoeffs.lookup(Cptr[11],0,t23,0,mq,mq,mq);
  PVcoeffs.lookup(Cptr[12],0,sqqb,s12,mq,mq,mq);


  PVcoeffs.lookup(Bptr[5], t13, mq,mq);
  PVcoeffs.lookup(Bptr[6], s12, mq,mq);
  PVcoeffs.lookup(Bptr[7], t23, mq,mq);
  PVcoeffs.lookup(Bptr[8], sqqb, mq,mq);


  
  /*
   * Remaining topologies
   */


  PVcoeffs.lookup(Bptr[9],0,0,0); 

  PVcoeffs.lookup(Bptr[18], t1q, 0,mq); 
  PVcoeffs.lookup(Bptr[19], t1qb, 0,mq); 
  PVcoeffs.lookup(Bptr[20], t2q, 0,mq); 
  PVcoeffs.lookup(Bptr[21], t2qb, 0,mq); 
  PVcoeffs.lookup(Bptr[22], s3q, 0,mq); 
  PVcoeffs.lookup(Bptr[23], s3qb,0,mq); 

  PVcoeffs.lookup(Bptr[24],t1qb, mq,0); 
  PVcoeffs.lookup(Bptr[25],t1q, mq,0); 
  PVcoeffs.lookup(Bptr[26], s3q, mq,0); 
  PVcoeffs.lookup(Bptr[27], s3qb, mq,0); 
  PVcoeffs.lookup(Bptr[28], t2q, mq,0); 
  PVcoeffs.lookup(Bptr[29], t2qb, mq,0); 


  /* 
   * sed s/"id C(C00\([ -z]*\)Cptr\([0-9]*\)(C00);"/"PVcoeffs.lookup(Cptr\2\1"/g tescht
   + Aufruf von Formskript...
   */

  PVcoeffs.lookup(Cptr[13],mq,t23,t1q,mq,0,0);
  PVcoeffs.lookup(Cptr[14],mq,s12,s3q,mq,0,0);
  PVcoeffs.lookup(Cptr[15],mq,t13,t2q,mq,0,0);
  PVcoeffs.lookup(Cptr[16],mq,t23,t1qb,mq,0,0);
  PVcoeffs.lookup(Cptr[17],mq,s12,s3qb,mq,0,0);
  PVcoeffs.lookup(Cptr[18],mq,t13,t2qb,mq,0,0);
  PVcoeffs.lookup(Cptr[19],t1q,t23,mq,mq,0,0);
  PVcoeffs.lookup(Cptr[20],t2qb,t13,mq,mq,0,0);
  PVcoeffs.lookup(Cptr[21],s3q,s12,mq,mq,0,0);
  PVcoeffs.lookup(Cptr[22],s3qb,s12,mq,mq,0,0);
  PVcoeffs.lookup(Cptr[23],t1qb,t23,mq,mq,0,0);
  PVcoeffs.lookup(Cptr[24],t2q,t13,mq,mq,0,0);
  PVcoeffs.lookup(Cptr[25],mq,sqqb,mq,mq,0,0);
  PVcoeffs.lookup(Cptr[26],t1qb,0,s3q,mq,0,0);
  PVcoeffs.lookup(Cptr[27],t1qb,0,t2q,mq,0,0);
  PVcoeffs.lookup(Cptr[28],t1q,0,s3qb,mq,0,0);
  PVcoeffs.lookup(Cptr[29],t1q,0,t2qb,mq,0,0);
  PVcoeffs.lookup(Cptr[30],t2qb,0,t1q,mq,0,0);
  PVcoeffs.lookup(Cptr[31],t2q,0,t1qb,mq,0,0);
  PVcoeffs.lookup(Cptr[32],mq,t1qb,0,0,mq,0);
  PVcoeffs.lookup(Cptr[33],mq,s3q,s12,0,mq,0);
  PVcoeffs.lookup(Cptr[34],mq,t2q,t13,0,mq,0);
  PVcoeffs.lookup(Cptr[35],mq,t2qb,0,0,mq,0);
  PVcoeffs.lookup(Cptr[36],mq,s3qb,0,0,mq,0);
  PVcoeffs.lookup(Cptr[37],mq,t1q,0,0,mq,0);
  PVcoeffs.lookup(Cptr[38],mq,s3qb,s12,0,mq,0);
  PVcoeffs.lookup(Cptr[39],mq,t2qb,t13,0,mq,0);
  PVcoeffs.lookup(Cptr[40],mq,t2q,0,0,mq,0);
  PVcoeffs.lookup(Cptr[41],mq,s3q,0,0,mq,0);
  PVcoeffs.lookup(Cptr[42],t1qb,s3q,0,0,mq,0);
  PVcoeffs.lookup(Cptr[43],t1qb,t2q,0,0,mq,0);
  PVcoeffs.lookup(Cptr[44],t1q,s3qb,0,0,mq,0);
  PVcoeffs.lookup(Cptr[45],t1q,t2qb,0,0,mq,0);
  PVcoeffs.lookup(Cptr[46],t2qb,s3q,0,0,mq,0);
  PVcoeffs.lookup(Cptr[47],t2q,s3qb,0,0,mq,0);
  PVcoeffs.lookup(Cptr[48],mq,t23,t1q,0,mq,mq);
  PVcoeffs.lookup(Cptr[49],mq,sqqb,mq,0,mq,mq);
  PVcoeffs.lookup(Cptr[50],mq,t23,t1qb,0,mq,mq);
  PVcoeffs.lookup(Cptr[51],s12,mq,s3q,mq,mq,0);
  PVcoeffs.lookup(Cptr[52],s12,mq,s3qb,mq,mq,0);
  PVcoeffs.lookup(Cptr[53],t23,mq,t1q,mq,mq,0);
  PVcoeffs.lookup(Cptr[54],t23,mq,t1qb,mq,mq,0);
  PVcoeffs.lookup(Cptr[55],t13,mq,t2qb,mq,mq,0);
  PVcoeffs.lookup(Cptr[56],t13,mq,t2q,mq,mq,0);
  PVcoeffs.lookup(Cptr[57],t13,mq,t2q,mq,mq,0);
  PVcoeffs.lookup(Cptr[58],0,t1qb,mq,mq,mq,0);
  PVcoeffs.lookup(Cptr[59],0,t1q,mq,mq,mq,0);
  PVcoeffs.lookup(Cptr[60],0,s3q,t2qb,mq,mq,0);
  PVcoeffs.lookup(Cptr[61],0,s3qb,t2q,mq,mq,0);
  PVcoeffs.lookup(Cptr[62],s12,s3q,mq,mq,mq,0);
  PVcoeffs.lookup(Cptr[63],s12,s3qb,mq,mq,mq,0);
  PVcoeffs.lookup(Cptr[64],0,t2q,t1qb,mq,mq,0);
  PVcoeffs.lookup(Cptr[65],0,t2qb,t1q,mq,mq,0);
  PVcoeffs.lookup(Cptr[66],0,t1q,t2qb,mq,mq,0);
  PVcoeffs.lookup(Cptr[67],0,t1qb,t2q,mq,mq,0);
  PVcoeffs.lookup(Cptr[68],0,s3qb,mq,mq,mq,0);
  PVcoeffs.lookup(Cptr[69],0,s3q,mq,mq,mq,0);
  PVcoeffs.lookup(Cptr[70],0,s3q,t1qb,mq,mq,0);
  PVcoeffs.lookup(Cptr[71],0,s3qb,t1q,mq,mq,0);
  PVcoeffs.lookup(Cptr[72],0,t2qb,mq,mq,mq,0);
  PVcoeffs.lookup(Cptr[73],0,t2q,mq,mq,mq,0);
  PVcoeffs.lookup(Cptr[74],t13,t2q,mq,mq,mq,0);
  PVcoeffs.lookup(Cptr[75],t13,t2qb,mq,mq,mq,0);
  PVcoeffs.lookup(Cptr[76],t23,t1q,mq,mq,mq,0);
  PVcoeffs.lookup(Cptr[77],t23,t1qb,mq,mq,mq,0);
  PVcoeffs.lookup(Cptr[78],0,t23,sqqb,0,0,0);
  PVcoeffs.lookup(Cptr[79],s12,sqqb,0,0,0,0);
  PVcoeffs.lookup(Cptr[80],s12,0,sqqb,0,0,0);
  PVcoeffs.lookup(Cptr[81],t13,0,sqqb,0,0,0);
  PVcoeffs.lookup(Cptr[82],0,s12,sqqb,0,0,0);
  PVcoeffs.lookup(Cptr[83],0,t13,sqqb,0,0,0);
  PVcoeffs.lookup(Cptr[84],s12,sqqb,0,mq,mq,mq);
  PVcoeffs.lookup(Cptr[85],s12,0,sqqb,mq,mq,mq);
  PVcoeffs.lookup(Cptr[86],t13,0,sqqb,mq,mq,mq);
  

  /*
   * sed s/"PVcoeffs.lookup(D001\([ -z]*\)Dptr\([0-9]*\)(D001);"/"PVcoeffs.lookup(Dptr\2\1"/g
   */

  PVcoeffs.lookup(Dptr[7],mq,t1qb,t23,sqqb,0,mq,0,mq,0,0);
  PVcoeffs.lookup(Dptr[8],mq,t2qb,t13,sqqb,0,mq,0,mq,0,0);
  PVcoeffs.lookup(Dptr[9],mq,s3qb,s12,sqqb,0,mq,0,mq,0,0);
  PVcoeffs.lookup(Dptr[10],mq,t1q,t23,sqqb,0,mq,0,mq,0,0);
  PVcoeffs.lookup(Dptr[11],mq,t2q,t13,sqqb,0,mq,0,mq,0,0);
  PVcoeffs.lookup(Dptr[12],mq,s3q,s12,sqqb,0,mq,0,mq,0,0);
  PVcoeffs.lookup(Dptr[13],t1qb,mq,s12,0,0,s3q,0,mq,0,0);
  PVcoeffs.lookup(Dptr[14],t1qb,mq,t13,0,0,t2q,0,mq,0,0);
  PVcoeffs.lookup(Dptr[15],t1q,mq,s12,0,0,s3qb,0,mq,0,0);
  PVcoeffs.lookup(Dptr[16],t1q,mq,t13,0,0,t2qb,0,mq,0,0);
  PVcoeffs.lookup(Dptr[17],t2qb,mq,t23,0,0,t1q,0,mq,0,0);
  PVcoeffs.lookup(Dptr[18],t2qb,s3q,s12,0,0,mq,0,mq,0,0);
  PVcoeffs.lookup(Dptr[19],t2q,mq,t23,0,0,t1qb,0,mq,0,0);
  PVcoeffs.lookup(Dptr[20],t2q,s3qb,s12,0,0,mq,0,mq,0,0);
  PVcoeffs.lookup(Dptr[21],s3qb,t2q,t13,0,0,mq,0,mq,0,0);
  PVcoeffs.lookup(Dptr[22],s3qb,t1q,t23,0,0,mq,0,mq,0,0);
  PVcoeffs.lookup(Dptr[23],s3q,t2qb,t13,0,0,mq,0,mq,0,0);
  PVcoeffs.lookup(Dptr[24],s3q,t1qb,t23,0,0,mq,0,mq,0,0);

  PVcoeffs.lookup(Dptr[25],mq,0,s3q,0,t2qb,t1qb,0,mq,mq,0);
  PVcoeffs.lookup(Dptr[26],mq,0,t2q,0,s3qb,t1qb,0,mq,mq,0);
  PVcoeffs.lookup(Dptr[27],mq,0,t1q,0,s3qb,t2qb,0,mq,mq,0);
  PVcoeffs.lookup(Dptr[28],mq,t23,mq,0,t1qb,t1q,0,mq,mq,0);
  PVcoeffs.lookup(Dptr[29],mq,t13,mq,0,t2qb,t2q,0,mq,mq,0);
  PVcoeffs.lookup(Dptr[30],mq,0,s3qb,0,t2q,t1q,0,mq,mq,0);
  PVcoeffs.lookup(Dptr[31],mq,s12,mq,0,s3qb,s3q,0,mq,mq,0);
  PVcoeffs.lookup(Dptr[32],mq,0,t2qb,0,s3q,t1q,0,mq,mq,0);
  PVcoeffs.lookup(Dptr[33],mq,0,t1qb,0,s3q,t2q,0,mq,mq,0);
  PVcoeffs.lookup(Dptr[34],0,t1qb,0,t2qb,mq,s3q,mq,mq,0,0);
  PVcoeffs.lookup(Dptr[35],0,t1qb,0,s3qb,mq,t2q,mq,mq,0,0);
  PVcoeffs.lookup(Dptr[36],0,t1q,t23,t1qb,mq,mq,mq,mq,0,0);
  PVcoeffs.lookup(Dptr[37],0,t1q,0,t2q,mq,s3qb,mq,mq,0,0);
  PVcoeffs.lookup(Dptr[38],0,t1q,0,s3q,mq,t2qb,mq,mq,0,0);
  PVcoeffs.lookup(Dptr[39],0,t2qb,0,s3qb,mq,t1q,mq,mq,0,0);
  PVcoeffs.lookup(Dptr[40],0,t2q,t13,t2qb,mq,mq,mq,mq,0,0);
  PVcoeffs.lookup(Dptr[41],0,t2q,0,s3q,mq,t1qb,mq,mq,0,0);
  PVcoeffs.lookup(Dptr[42],0,s3q,s12,s3qb,mq,mq,mq,mq,0,0);
  
  PVcoeffs.lookup(Dptr[43],0,t23,mq,mq,sqqb,t1q,mq,mq,mq,0);
  PVcoeffs.lookup(Dptr[44],0,t23,mq,mq,sqqb,t1qb,mq,mq,mq,0);
  PVcoeffs.lookup(Dptr[45],0,s12,mq,t2qb,0,s3q,mq,mq,mq,0);
  PVcoeffs.lookup(Dptr[46],0,s12,mq,t2q,0,s3qb,mq,mq,mq,0);
  PVcoeffs.lookup(Dptr[47],0,s12,s3q,t1qb,0,mq,mq,mq,mq,0);
  PVcoeffs.lookup(Dptr[48],0,s12,s3qb,t1q,0,mq,mq,mq,mq,0);
  PVcoeffs.lookup(Dptr[49],0,t13,mq,s3qb,0,t2q,mq,mq,mq,0);
  PVcoeffs.lookup(Dptr[50],0,t13,mq,s3q,0,t2qb,mq,mq,mq,0);
  PVcoeffs.lookup(Dptr[51],0,t13,t2q,t1qb,0,mq,mq,mq,mq,0);
  PVcoeffs.lookup(Dptr[52],0,t13,t2qb,t1q,0,mq,mq,mq,mq,0);
  PVcoeffs.lookup(Dptr[53],0,t13,mq,mq,sqqb,t2q,mq,mq,mq,0);
  PVcoeffs.lookup(Dptr[54],0,t13,mq,mq,sqqb,t2qb,mq,mq,mq,0);
  PVcoeffs.lookup(Dptr[55],0,t23,mq,s3qb,0,t1q,mq,mq,mq,0);
  PVcoeffs.lookup(Dptr[56],0,t23,mq,s3q,0,t1qb,mq,mq,mq,0);
  PVcoeffs.lookup(Dptr[57],0,t23,t1q,t2qb,0,mq,mq,mq,mq,0);
  PVcoeffs.lookup(Dptr[58],0,t23,t1qb,t2q,0,mq,mq,mq,mq,0);
  PVcoeffs.lookup(Dptr[59],0,s12,mq,mq,sqqb,s3q,mq,mq,mq,0);
  PVcoeffs.lookup(Dptr[60],0,s12,mq,mq,sqqb,s3qb,mq,mq,mq,0);

  
}


double gggtt2virt(const FourMomentum & p1,const FourMomentum & p2,
		  const FourMomentum & p3,
		  const FourMomentum & kq, const FourMomentum & kqb){

  static ErrorCounter LargeCorrection("Large correction drop event.");
  FourMomentum q1, q2, r1, r2; 

  const double N = 3.0;
  const double CF = 1.0/2.0/N*(N*N-1.0);
  const double nfl = 5.0;
  const double m = kq.getmass();
  const double mq = kq.getmass2();


  complex<double> ampLO[7] ;
  double amp2 = 0;
  double amplo2 = 0;
  complex<double> ampNLO[12] ;
  complex<double>S[2][2][2][133] ;
  complex<double> Mij[7][7] ;
  
  IntType f[12][133] ;
  double flo[7][133] ;
  double mct[7][133];
  IntType f_nfl[12][133];
  IntType f_top[12][133];
  
  IntType fivepoint[12][133];
  IntType rest[12][133];
  IntType ghosts[12][133];
  
  complex<double> tmp;
  
  KahanSum summer;




  static bool initialized = false;

  static double CNLO[12][12] __attribute__((aligned(64)));

  CoeffCache PVcoeffs;
  IntType* Bptr[30]={0};
  IntType* Cptr[87]={0};
  IntType* Dptr[61]={0};


  if ( !initialized ) {
    /*
     * Evaluate the color correlation matrix
     */ 
    cout << "# -------------------------------------------------------------"
	 << endl;
    cout << "# gggtt2virt: Initialize color matrix..."; 
    NLOcolorgg(CNLO);
    initialized = true;
    cout << " done." << endl;
    cout << "# -------------------------------------------------------------"
	 << endl;
  }

  /*
   * Start with the evaluation of the pentagons, 
   * if an exception is thrown we leave the function:
   */
  try {
    pentagon_gg(fivepoint,p1,p2,p3,kq,kqb, 0.0); 
  } catch (GGException & ex) {
    ex.print();
    Topology::cacheclear();
    return(0.);
  } catch ( logic_error & ex) {
    cout << "# " << ex.what() << endl;
    Topology::cacheclear();
    return(0.);
  }

  evaluatePVcoefficients_gg(PVcoeffs,p1,p2,p3,kq,kqb,Bptr,Cptr,Dptr);

  MasslessFermionLoops_gg(f_nfl,p1,p2,p3,kq,kqb,Bptr,Cptr,Dptr,1.0 );   

  TopLoops_gg(f_top,p1,p2,p3,kq,kqb,Bptr,Cptr,Dptr, 1.0); 

  RemainingDiags_gg(rest,p1,p2,p3,kq,kqb,Bptr,Cptr,Dptr, 1.0);

  Ghosts_gg(ghosts,p1,p2,p3,kq,kqb,Bptr,Cptr,Dptr, 1.0);

  masscounterterm_gg(mct,p1,p2,p3,kq,kqb, 1.0);

  EvalLOfgg(flo, p1, p2, p3, kq, kqb);
 
  /*
   * Sum up the different contibutions to NLO corrections
   */

  for( int color=1; color < 12; color++ ){
    for( int i=1; i < 133; i++ ){
      f[color][i] =
	+ nfl * f_nfl[color][i]
	+ f_top[color][i]
	+ rest[color][i] 
	+ ghosts[color][i] 
	+ fivepoint[color][i]
	;
    }
  }

  iopgg(Mij,p1,p2,p3,kq,kqb);

  /*
   * Add the contribution from the renormalization:
   */
  const double RenormalizationScale2 = 
    parms.getRenormalizationScaleSquared(); 
  const double InternalScale2 = parms.getInternalScaleSquared();
  const double deltaUV1 = parms.getDeltaUV1() + log(InternalScale2/mq);
  const double deltaIR1 = parms.getDeltaIR1() + log(InternalScale2/mq);

  /*
   * Renormalization constants, factor alpha_s/4/pi is taken out
   */
  const double dm_over_m = - CF * ( 3.0*deltaUV1 + 4.0 );

  const double dalphas_over_alphas = 
    ( (2.0/3.0) * nfl - (11.0/3.0) * N ) 
    * (deltaUV1 + log(mq/RenormalizationScale2))
    + (2.0/3.0) * deltaUV1;

  const double dZt = - CF * ( deltaUV1 + 2.0 * deltaIR1 + 4.0 );

  const double dZA = 
    - ( (2.0/3.0) * nfl   - (5.0/3.0) * N ) * ( deltaUV1 - deltaIR1 )
    - (2.0/3.0) * deltaUV1;

  const double born_renormalization =  
    (3.0/2.0) * dZA + dZt + (3.0/2.0) * dalphas_over_alphas ;

  for(int color=1; color < 7; color++){
    for(int i=1; i < 133; i++){
       f[color][i] += 
	born_renormalization * flo[color][i]
	+ dm_over_m * mct[color][i];
    }
  }


  evalqr(kq,kqb,q1,q2,r1,r2);

  /*
   *  Sum over the polarization and square it:
   */

  EvalSMEgg(S,p1, p2, p3,q1,q2,r1,r2);

  for (int tspin=0; tspin < 2; tspin++) {

    if (tspin == 1) 
      EvalSMEgg(S,p1, p2, p3,q1,q2,r2,r1);

    //if (tspin == 2) 
    //  EvalSMEgg(S,p1, p2, p3,q2,q1,r1,r2);

    //    if (tspin == 3) 
    //  EvalSMEgg(S,p1, p2, p3,q2,q1,r2,r1);

    for (int pol1=0; pol1 < 2; pol1++){
      for (int pol2=0; pol2 < 2; pol2++){
	for (int pol3=0; pol3 < 2; pol3++){
	  /*
	   *  Calculate the leading order amplitude, the
           *  result is vector a in color space 
	   */
	  for(int color=1; color < 7; color++){
	    summer.reset();
	    for(int i=1;i<133;i++){
	      //   ampLO[color] += S[pol1][pol2][pol3][i] * flo[color][i];
	      summer.add( S[pol1][pol2][pol3][i] * flo[color][i] );
	    }
	    ampLO[color] = summer.getSum();
	  }

	  /*
	   *  Calculate the next-to-leading order amplitude, the
           *  result is a vector in color space 
	   */
	  for(int color=1; color < 12; color++){
	    summer.reset();
	    for(int i=1; i < 133; i++){
	      // ampNLO[color] += S[pol1][pol2][pol3][i] * f_nfl[color][i];
	      summer.add(S[pol1][pol2][pol3][i] * f[color][i] );
	    }
	    ampNLO[color] = summer.getSum();
	  }

	  /*
	   *  Calculate the inference term between leading oder
	   *  and next-to-leading order amplitude
	   */
	  summer.reset();
	  for(int c1=1; c1 < 7; c1++) {
	    for(int c2=1; c2 < 12; c2++) {
	      //tmp += conj(ampLO[c1]) * CNLO[c1][c2] * ampNLO[c2];
	      summer.add( conj(ampLO[c1]) * CNLO[c1][c2] * ampNLO[c2] ); 
	    }
	  }

	  /*
	   * Add the contribution from the Ioperator:
	   */
	  complex<double> borntmp(0.,0.);
	  for(int c1=1; c1 < 7; c1++) {
	    for(int c2=1; c2 < 7; c2++) {
	      summer.add( conj(ampLO[c1]) * Mij[c1][c2] * ampLO[c2] ); 
	      borntmp +=  conj(ampLO[c1]) * CNLO[c1][c2] * ampLO[c2];
	    }
	  }

	  amp2 += summer.getSum().real();
	  amplo2 += borntmp.real();
	} // pol3 loop
      } // pol2 loop
    } // pol1 loop
  } // tspin loop
#ifdef TESTBADPOINTS
  static EPSHisto epsgg("EPS_large_gg.top");
  epsgg.fill(log(fabs(0.01*amp2/amplo2)));
#endif
  if ( fabs(0.01*amp2/amplo2) > 100. ) {
    LargeCorrection.addError();
    amp2 = 0.0;
  }
  // 2 (spin) * 2 (2Re()):
  return (4.*amp2);
}

