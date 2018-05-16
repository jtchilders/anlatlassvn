// $Modified: Tue Jul 31 15:29:22 2007 by uwer $
#include "StandardModelParameters.h"
#include "EvalSMExx.h"
#include "FourMomentum.h"
#include "Subtractions.h"
#include "Color.h"
#include "Reduction.h"
#include "GGReduction.h"
#include "ScalarInt.h"
#include "Evalf.h"
#include <complex>
#include <cmath>
#include <cstring>
#include <iostream>
#include "EPSHisto.h"

using namespace std;

static StandardModelParameters& parms = StandardModelParameters::instance();


#define MAXCOLOR 5


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


void iopqg(complex<double> Mij[MAXCOLOR][MAXCOLOR],
		    const FourMomentum & p1, const FourMomentum & p2, 
		    const FourMomentum & p3,
		    const FourMomentum & kt, const FourMomentum & ktb ) {

  double m = kt.getmass();

  static double clo[5][5];
  static double cqqb[5][5];
  static double cttb[5][5];
  static double cqg3[5][5]; 
  static double cqbg3[5][5]; 
  static double cqt[5][5];
  static double cqbt[5][5]; 
  static double cg3t[5][5]; 
  static double cqtb[5][5];
  static double cqbtb[5][5]; 
  static double cg3tb[5][5];

  static bool initialized = false;

  double mu2 = parms.getRenormalizationScaleSquared();

  if (!initialized) {
    qqcorrelations(clo, cqqb, cttb,cqg3, cqbg3, cqt,
		   cqbt, cg3t, cqtb,cqbtb, cg3tb);
    initialized = true;
  }


  const double kappa = (2.0/3.0);
  vector<double> mf;

  for(int c1=1; c1 < MAXCOLOR; c1++) {
    for(int c2=c1; c2 < MAXCOLOR; c2++) {

      using namespace Subtractions;

      /*
       * In comparison with the ggttg process the main difference comes
       * from the "initial state emitter, final state spectator" configuration
       * Apart from that only the color factors are changed.
       */

      Mij[c1][c2] =  -1.0 
	* (
	   // Final state I_m
	   + cttb[c1][c2] 
	   * ( Vq( 2.0*dotp(kt,ktb),m,m) - (pi2/3.0) 
	       + ( Gamma_q(m) + gammaq*log(mu2/2.0/dotp(kt,ktb)) 
		   + gammaq + Kq )/CF 
	       )
	   
	   + cttb[c1][c2] 
	   * ( Vq( 2.0*dotp(kt,ktb),m,m) - (pi2/3.0) 
	       + ( Gamma_q(m) + gammaq*log(mu2/2.0/dotp(kt,ktb)) 
		   + gammaq + Kq )/CF 
	       )
	   
	   + cqbt[c1][c2] 
	   * ( Vq( 2.0*dotp(kt,p3),m,0.) - (pi2/3.0) 
	       + ( Gamma_q(m) + gammaq*log(mu2/2.0/dotp(kt,p3)) 
		   + gammaq + Kq)/CF 
	       ) 
	   
	   + cqbtb[c1][c2] 
	   * ( Vq( 2.0*dotp(ktb,p3),m,0.) - (pi2/3.0) 
	       + ( Gamma_q(m) + gammaq*log(mu2/2.0/dotp(ktb,p3)) 
		   + gammaq + Kq)/CF 
	       )
	   
	   + cqbt[c1][c2] 
	   * ( Vq( 2.0*dotp(p3,kt), 0., m) - (pi2/3.0) 
	       + ( Gamma_q() + gammaq*log(mu2/2.0/dotp(p3,kt)) 
		   + gammaq + Kq )/CF 
	       )
	   
	   + cqbtb[c1][c2] 
	   * ( Vq( 2.0*dotp(p3,ktb), 0., m) - (pi2/3.0) 
	       + ( Gamma_q() + gammaq*log(mu2/2.0/dotp(p3,ktb)) 
		   + gammaq + Kq )/CF 
	       )
	   
	   // initial -- final I_a 
	   /*
	    * Final state emitter, initial state spectator
	    */
	   + cqt[c1][c2] 
	   * ( Vq( 2.0*dotp(kt,p1),m,0.) - (pi2/3.0) 
	       + ( Gamma_q(m) + gammaq*log(mu2/2.0/dotp(kt,p1)) 
		   + gammaq + Kq )/CF 
	       )
	   
	   + cqtb[c1][c2] 
	   * ( Vq( 2.0*dotp(ktb,p1),m,0.) - (pi2/3.0) 
	       + ( Gamma_q(m) + gammaq*log(mu2/2.0/dotp(ktb,p1)) 
		   + gammaq + Kq )/CF 
	       )
	   
	   + cqqb[c1][c2] 
	   * ( Vq( 2.0*dotp(p3,p1), 0., 0.) - (pi2/3.0) 
	       + ( Gamma_q() + gammaq*log(mu2/2.0/dotp(p3,p1)) 
		   + gammaq + Kq )/CF 
	       )
	   /*
	    * Final state spectator, initial state emitter
	    */
	   + cqt[c1][c2] 
	   * ( Vq( 2.0*dotp(p1,kt), 0., m) - (pi2/3.0) 
	       + ( Gamma_q() + gammaq*log(mu2/2.0/dotp(p1,kt)) 
		   + gammaq + Kq )/CF 
	       )
	   
	   + cqtb[c1][c2] 
	   * ( Vq( 2.0*dotp(p1,ktb), 0., m) - (pi2/3.0) 
	       + ( Gamma_q() + gammaq*log(mu2/2.0/dotp(p1,ktb)) 
		   + gammaq + Kq )/CF 
	       )
	   + cqqb[c1][c2] 
	   * ( Vq( 2.0*dotp(p1,p3), 0., 0.) - (pi2/3.0) 
	       + ( Gamma_q() + gammaq*log(mu2/2.0/dotp(p1,p3)) 
		   + gammaq + Kq )/CF 
	       )
	   // initial -- final I_b 
	   /*
	    * Final state emitter, initial state spectator
	    */
	   + cg3t[c1][c2] 
	   * ( Vq( 2.0*dotp(kt,p2),m,0.) - (pi2/3.0) 
	       + ( Gamma_q(m) + gammaq*log(mu2/2.0/dotp(kt,p2)) 
		   + gammaq + Kq )/CF 
	       )
	   
	   + cg3tb[c1][c2] 
	   * ( Vq( 2.0*dotp(ktb,p2),m,0.) - (pi2/3.0) 
	       + ( Gamma_q(m) + gammaq*log(mu2/2.0/dotp(ktb,p2)) 
		   + gammaq + Kq)/CF 
	       )
	   
	   + cqbg3[c1][c2] 
	   * ( Vq( 2.0*dotp(p3,p2), 0., 0.) - (pi2/3.0) 
	       + ( Gamma_q() + gammaq*log(mu2/2.0/dotp(p3,p2)) 
		   + gammaq + Kq )/CF 
	       )
	   /*
	    * Final state spectator, initial state emitter
	    */
	   + cg3t[c1][c2] 
	   * ( Vg( 2.0*dotp(p2,kt), 0., m,mf, kappa) - (pi2/3.0) 
	   + ( Gamma_g(mf) + gammag*log(mu2/2.0/dotp(p2,kt)) 
	       + gammag + Kg )/CA 
	       )
	   
	   + cg3tb[c1][c2] 
	   * ( Vg( 2.0*dotp(p2,ktb), 0., m,mf, kappa) - (pi2/3.0) 
	       + ( Gamma_g(mf) + gammag*log(mu2/2.0/dotp(p2,ktb)) 
		   + gammag + Kg )/CA 
	       )
	   + cqbg3[c1][c2] 
	   * ( Vg( 2.0*dotp(p2,p3), 0., 0.,mf, kappa) - (pi2/3.0) 
	       + ( Gamma_g(mf) + gammag*log(mu2/2.0/dotp(p2,p3)) 
		   + gammag + Kg )/CA 
	       )
	   /*
	    * The two additional contributions from  
	    * [hep-ph/0201036, eq. 6.66] 
	    */ 
	   + cqg3[c1][c2] 
	   * ( Vsing(2.0*dotp(p1,p2), 0., 0.) - (pi2/3.0) 
	       + ( Gamma_q() + gammaq * log(mu2/2.0/dotp(p1,p2)) 
		   + gammaq + Kq )/CF
	       )
	   
	   + cqg3[c1][c2] 
	   * ( Vsing(2.0*dotp(p1,p2), 0., 0.) - (pi2/3.0) 
	       + ( Gamma_g() + gammag * log(mu2/2.0/dotp(p1,p2)) 
		   + gammag + Kg )/CA
	       )
	   );   
      Mij[c2][c1] = Mij[c1][c2];
    }
  }
}
    
double Ioperator_qg(const FourMomentum & p1, const FourMomentum & p2, 
		 const FourMomentum & p3,
		 const FourMomentum & kq, const FourMomentum & kqb){

  FourMomentum q1, q2, r1, r2; 


  complex<double>S[2][2][SMEMAXQQ];
  double f[MAXCOLOR][SMEMAXQQ];
  complex<double> ampLO[MAXCOLOR];

  complex<double> Mij[MAXCOLOR][MAXCOLOR],tmp;

  double amp2=0;


  static bool initialized=false;


  evalqr(kq,kqb,q1,q2,r1,r2);
  EvalLOfqq(f, p1,-p3,-p2,kq,kqb);

  iopqg(Mij,p1,p2,p3,kq,kqb);

  /*
   *  Sum over the polarization and square it:
   */
  EvalSMEqq(S,p1, -p3, -p2,q1,q2,r1,r2);

  for (int tspin=0; tspin < 2; tspin++) {

    if (tspin == 1) 
      EvalSMEqq(S,p1, -p3, -p2,q1,q2,r2,r1);

    if (tspin == 2) 
      EvalSMEqq(S,p1, -p3, -p2,q2,q1,r1,r2);

    if (tspin == 3) 
      EvalSMEqq(S,p1, -p3, -p2,q2,q1,r2,r1);

    for (int pol1=0; pol1 < 2; pol1++){
      for (int pol3=0; pol3 < 2; pol3++){
	for(int color=1;color < MAXCOLOR;color++){
	  ampLO[color] = 0.0;
	  for(int i=1; i < SMEMAXQQ; i++){
	    ampLO[color] +=  S[pol1][pol3][i] * f[color][i];
	  }
	}
	tmp = complex<double>(0.,0.);
	for(int c1=1; c1 < MAXCOLOR; c1++){
	  for(int c2=1; c2 < MAXCOLOR; c2++){
	    tmp += conj(ampLO[c1]) * Mij[c1][c2] * ampLO[c2];
	  }
	}
	amp2 = amp2 + tmp.real();
      }
    }
  }
  return(2.0*amp2);
}

double qgttq2(const FourMomentum & p1, const FourMomentum & p2, 
              const FourMomentum & p3,
              const FourMomentum & kq, const FourMomentum & kqb){

  /*
   * q(p1) g(p2) --> t(kq) tb(kqb) + q(p3) is obtained from
   * 
   * q(p1) qb(p2) --> t(kq) tb(kqb) + g(p3) by crossing
   *
   * p2 --> - p3, p3 --> -p2
   *
   */

  FourMomentum q1("q1"), q2("q2"), r1("r1"), r2("r2"); 

  complex<double> S[2][2][SMEMAXQQ];

  double f[MAXCOLOR][SMEMAXQQ];  

  complex<double> ampLO[MAXCOLOR], tmp;

  double amp2=0;

  static double CLO[MAXCOLOR][MAXCOLOR];

  static bool initialized=false;


  if ( !initialized ) {
    LOcolorqq(CLO);
    initialized = true;
  }

  evalqr(kq,kqb,q1,q2,r1,r2);
  EvalLOfqq(f, p1,-p3,-p2,kq,kqb);



  /*
   *  Sum over the polarization and square it:
   */

  EvalSMEqq(S,p1, -p3, -p2,q1,q2,r1,r2);

  for (int tspin=0; tspin < 2; tspin++) {

    if (tspin == 1) 
      EvalSMEqq(S,p1, -p3, -p2,q1,q2,r2,r1);    

    if (tspin == 2) 
      EvalSMEqq(S,p1, -p3, -p2,q2,q1,r1,r2);

    if (tspin == 3) 
      EvalSMEqq(S,p1, -p3, -p2,q2,q1,r2,r1);

    for (int pol1=0; pol1 < 2; pol1++){
      for (int pol3=0; pol3 < 2; pol3++){
	
	for(int color=1;color < MAXCOLOR;color++){
	  ampLO[color] = complex<double>(0.,0.);
	  for(int i=1; i < SMEMAXQQ; i++){
	    ampLO[color] +=  S[pol1][pol3][i] * f[color][i];
	  }
	}

	tmp = complex<double>(0.,0.);
	for(int c1=1; c1 < MAXCOLOR; c1++) {
	  for(int c2=1; c2 < MAXCOLOR; c2++) {
	    tmp = tmp + conj(ampLO[c1]) * CLO[c1][c2] * ampLO[c2];
	  }
	}
	amp2 = amp2 + tmp.real();
      }
    }
  }
  return 2.0 * amp2;
}


static void evaluatePVcoefficients(CoeffCache& PVcoeffs,
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


  PVcoeffs.lookup(Bptr[1],s12,0,0);
  PVcoeffs.lookup(Bptr[2],t13,0,0);
  PVcoeffs.lookup(Bptr[3],sqqb,0,0);
  PVcoeffs.lookup(Bptr[4],t23,0,0);
  PVcoeffs.lookup(Bptr[5],s3q,0,mq);
  PVcoeffs.lookup(Bptr[6],s3qb,0,mq);
  //  PVcoeffs.lookup(Bptr[7],sqqb,0,0);
  PVcoeffs.lookup(Bptr[8],0,0,0);
  PVcoeffs.lookup(Bptr[9],s12,mq,mq);
  PVcoeffs.lookup(Bptr[10],sqqb,mq,mq);
  // PVcoeffs.lookup(Bptr[11],sqqb,mq,mq);
  

  PVcoeffs.lookup(Cptr[1],0,sqqb,s12,mq,mq,mq);
  PVcoeffs.lookup(Cptr[2],s12,0,sqqb,0,0,0);
  PVcoeffs.lookup(Cptr[3],s12,sqqb,0,0,0,0);
  PVcoeffs.lookup(Cptr[4],s12,s3qb,mq,mq,mq,0);
  PVcoeffs.lookup(Cptr[5],s12,s3q,mq,mq,mq,0);
  PVcoeffs.lookup(Cptr[6],s3q,s12,mq,mq,0,0);
  PVcoeffs.lookup(Cptr[7],s3qb,s12,mq,mq,0,0);
  PVcoeffs.lookup(Cptr[8],0,sqqb,s12,0,0,0);
  PVcoeffs.lookup(Cptr[9],0,s3q,mq,mq,mq,0);
  PVcoeffs.lookup(Cptr[10],0,s3qb,mq,mq,mq,0);
  PVcoeffs.lookup(Cptr[11],mq,s3qb,s12,0,mq,0);
  PVcoeffs.lookup(Cptr[12],mq,sqqb,mq,mq,0,0);
  PVcoeffs.lookup(Cptr[13],mq,sqqb,mq,0,mq,mq);
  PVcoeffs.lookup(Cptr[14],mq,s3q,s12,0,mq,0);
  PVcoeffs.lookup(Cptr[15],0,s12,0,0,0,0);
  PVcoeffs.lookup(Cptr[16],0,sqqb,t23,0,0,0);
  PVcoeffs.lookup(Cptr[17],0,t13,0,0,0,0);
  PVcoeffs.lookup(Cptr[18],0,sqqb,t13,0,0,0);
  PVcoeffs.lookup(Cptr[19],0,t23,0,0,0,0);
  PVcoeffs.lookup(Cptr[20],0,t13,0,0,0,0);
  PVcoeffs.lookup(Cptr[21],0,t23,0,0,0,0);
  PVcoeffs.lookup(Cptr[22],mq,s3q,0,0,mq,0);
  PVcoeffs.lookup(Cptr[23],mq,s3qb,0,0,mq,0);



  PVcoeffs.lookup(Dptr[1],t1q,mq,s12,0,0,s3qb,0,mq,0,0);
  PVcoeffs.lookup(Dptr[2],t1qb,mq,s12,0,0,s3q,0,mq,0,0);
  PVcoeffs.lookup(Dptr[3],t2q,s3qb,s12,0,0,mq,0,mq,0,0);
  PVcoeffs.lookup(Dptr[4],t2qb,s3q,s12,0,0,mq,0,mq,0,0);
  PVcoeffs.lookup(Dptr[5],0,s12,0,t23,0,sqqb,0,0,0,0);
  PVcoeffs.lookup(Dptr[6],0,sqqb,0,0,t23,t13,0,0,0,0);
  PVcoeffs.lookup(Dptr[7],0,t13,sqqb,t23,0,0,0,0,0,0);
  PVcoeffs.lookup(Dptr[8],0,s3q,s12,s3qb,mq,mq,mq,mq,0,0);
  PVcoeffs.lookup(Dptr[9],mq,s3q,s12,sqqb,0,mq,0,mq,0,0);
  PVcoeffs.lookup(Dptr[10],mq,mq,t23,0,sqqb,t1q,0,mq,0,0);
  PVcoeffs.lookup(Dptr[11],mq,mq,t13,0,sqqb,t2q,0,mq,0,0);
  PVcoeffs.lookup(Dptr[12],mq,s3qb,s12,sqqb,0,mq,0,mq,0,0);
  PVcoeffs.lookup(Dptr[13],mq,mq,t23,0,sqqb,t1qb,0,mq,0,0);
  PVcoeffs.lookup(Dptr[14],mq,mq,t13,0,sqqb,t2qb,0,mq,0,0);
  PVcoeffs.lookup(Dptr[15],0,s12,sqqb,t13,0,0,0,0,0,0);
  PVcoeffs.lookup(Dptr[16],0,s12,mq,mq,sqqb,s3qb,mq,mq,mq,0);
  PVcoeffs.lookup(Dptr[17],0,s12,mq,mq,sqqb,s3q,mq,mq,mq,0);
  PVcoeffs.lookup(Dptr[18],mq,s12,mq,0,s3qb,s3q,0,mq,mq,0);


}


double qgttq2virt(const FourMomentum & p1,const FourMomentum & p2,
		  const FourMomentum & p3,
		  const FourMomentum & kq, const FourMomentum & kqb){


  FourMomentum q1, q2, r1, r2; 

  const double N = 3.0;
  const double CF = (4.0/3.0);
  const double nfl = 5.0;
  const double m = kq.getmass();
  const double mq = kq.getmass2();

  // __attribute__((aligned(64)))
  complex<double> ampLO[MAXCOLOR] ;
  double amp2 = 0;
  double amplo2 = 0;
  complex<double> ampNLO[MAXCOLOR] ;
  complex<double>S[2][2][SMEMAXQQ] ;
  IntType f[MAXCOLOR][SMEMAXQQ] ;
  double flo[MAXCOLOR][SMEMAXQQ] ;
  double mct[MAXCOLOR][SMEMAXQQ];
  IntType f_nfl[MAXCOLOR][SMEMAXQQ];
  IntType f_top[MAXCOLOR][SMEMAXQQ];

  IntType fivepoint[MAXCOLOR][SMEMAXQQ];
  IntType rest[MAXCOLOR][SMEMAXQQ];
  IntType ghosts[MAXCOLOR][SMEMAXQQ];

  complex<double> Mij[MAXCOLOR][MAXCOLOR];

  KahanSum summer;

  complex<double> tmp;


  static bool initialized = false;

  static double CLO[MAXCOLOR][MAXCOLOR];

  CoeffCache PVcoeffs;
  IntType* Bptr[12]={0};
  IntType* Cptr[24]={0};
  IntType* Dptr[19]={0};


  if ( !initialized ) {
    /*
     * Evaluate the color correlation matrix
     */ 
    cout << "# -------------------------------------------------------------"
	 << endl;
    cout << "# qgttq2virt: Initialize color matrix..."; 
    LOcolorqq(CLO);
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
    pentagon_qq(fivepoint,p1,-p3,-p2,kq,kqb, 0.0); 
  } catch (GGException & ex) {
    ex.print();
    Topology::cacheclear();
    return(0.);
  } catch ( logic_error & ex) {
    cout << "# " << ex.what() << endl;
    Topology::cacheclear();
    return(0.);
  }


  evaluatePVcoefficients(PVcoeffs,p1,-p3,-p2,kq,kqb,Bptr,Cptr,Dptr);

  MasslessFermionLoops_qq(f_nfl,p1,-p3,-p2,kq,kqb,Bptr,Cptr,Dptr,1.0 );   

  TopLoops_qq(f_top,p1,-p3,-p2,kq,kqb,Bptr,Cptr,Dptr, 1.0); 

  RemainingDiags_qq(rest,p1,-p3,-p2,kq,kqb,Bptr,Cptr,Dptr, 1.0);

  Ghosts_qq(ghosts,p1,-p3,-p2,kq,kqb,Bptr,Cptr,Dptr, 1.0);

  masscounterterm_qq(mct,p1,-p3,-p2,kq,kqb, 1.0);

  EvalLOfqq(flo, p1, -p3, -p2, kq, kqb);

  iopqg(Mij,p1,p2,p3,kq,kqb);
  /*
   * Sum up the different contibutions to NLO corrections
   */

  for( int color=1; color < MAXCOLOR; color++ ){
    for( int i=1; i < SMEMAXQQ; i++ ){
      f[color][i] = 0.0
	+ nfl * f_nfl[color][i]
	+ f_top[color][i]
	+ rest[color][i] 
	+ ghosts[color][i] 
	+ fivepoint[color][i]
	;
    }
  }

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
  const double fermionic = 1.0;
  const double bosonic = 1.0;

  const double dm_over_m = - CF * ( 3.0*deltaUV1 + 4.0 ) * bosonic;

  const double dalphas_over_alphas = 
    ( (2.0/3.0) * nfl * fermionic - (11.0/3.0) * N * bosonic ) 
    * ( deltaUV1 + log(mq/RenormalizationScale2) )
    + (2.0/3.0) * deltaUV1 * fermionic;

  const double dZt = - CF * ( deltaUV1 + 2.0 * deltaIR1 + 4.0 ) * bosonic ;

  const double dZq = - CF * ( deltaUV1 - deltaIR1 ) * bosonic;

  const double dZA = 
    - ( + (2.0/3.0) * nfl * fermionic  
	- (5.0/3.0) * N  * bosonic ) * ( deltaUV1 - deltaIR1 )
    - (2.0/3.0) * deltaUV1 * fermionic;

  const double born_renormalization = 
    (1.0/2.0)*dZA + dZt + dZq + (3.0/2.0) * dalphas_over_alphas;


  for(int color=1; color < MAXCOLOR; color++){
    for(int i=1; i < SMEMAXQQ; i++){
      f[color][i] += 
	born_renormalization * flo[color][i]
	+ dm_over_m * mct[color][i];
    }
  }

  evalqr(kq,kqb,q1,q2,r1,r2);

  /*
   *  Sum over the polarization and square it:
   */


  EvalSMEqq(S,p1, -p3, -p2,q1,q2,r1,r2);

  for (int tspin=0; tspin < 2; tspin++) {

    if (tspin == 1) 
      EvalSMEqq(S,p1, -p3, -p2,q1,q2,r2,r1);

    //if (tspin == 2) 
    //  EvalSMEqq(S,p1, -p3, -p2,q2,q1,r1,r2);

    //    if (tspin == 3) 
    //  EvalSMEqq(S,p1, -p3, -p2,q2,q1,r2,r1);

    for (int pol1=0; pol1 < 2; pol1++){
	for (int pol3=0; pol3 < 2; pol3++){
	  /*
	   *  Calculate the leading order amplitude, the
           *  result is vector a in color space 
	   */
	  for(int color=1; color < MAXCOLOR; color++){
	    summer.reset();
	    for(int i=1;i<SMEMAXQQ;i++){
	      //   ampLO[color] += S[pol1][pol2][pol3][i] * flo[color][i];
	      summer.add( S[pol1][pol3][i] * flo[color][i] );
	    }
	    ampLO[color] = summer.getSum();
	  }

	  /*
	   *  Calculate the next-to-leading order amplitude, the
           *  result is a vector in color space 
	   */
	  for (int color=1; color < MAXCOLOR; color++){
	    summer.reset();
	    for(int i=1; i < SMEMAXQQ; i++){
	      // ampNLO[color] += S[pol1][pol2][pol3][i] * f_nfl[color][i];
	      summer.add(S[pol1][pol3][i] * f[color][i] );
	    }
	    ampNLO[color] = summer.getSum();
	  }

	  //	  cout << (ampNLO[1]+ampNLO[2]+ampNLO[3])/ampNLO[4] << endl;
	  /*
	   *  Calculate the inference term between leading oder
	   *  and next-to-leading order amplitude
	   *  plus contribution from Ioperator
	   */
	  summer.reset();
	  complex<double> borntmp(0.,0.);
	  for(int c1=1; c1 < MAXCOLOR; c1++) {
	    for(int c2=1; c2 < MAXCOLOR; c2++) {
	      //tmp += conj(ampLO[c1]) * CLO[c1][c2] * ampNLO[c2];
	      summer.add( conj(ampLO[c1]) * CLO[c1][c2] * ampNLO[c2] ); 
	      summer.add( conj(ampLO[c1]) * Mij[c1][c2] * ampLO[c2] ); 
	      borntmp += conj(ampLO[c1]) * CLO[c1][c2] * ampLO[c2];
	    }
	  }
	  amp2 += summer.getSum().real();
	  amplo2 += borntmp.real();
	} // pol3 loop
    } // pol1 loop
  } // tspin loop
#ifdef TESTBADPOINTS
  static EPSHisto epsqg("EPS_large_qg.top");
  epsqg.fill(log(fabs(0.01*amp2/amplo2)));
#endif
  if ( fabs(0.01*amp2/amplo2) > 100. ) {
    cout << "# Large correction drop event \n";
    cout << "# born = " << amplo2 << " nlo_qg = " << amp2*0.01 << endl;
    amp2 = 0.0;
  }
  // 2 (spin) * 2 (2Re())
  return(4.0*amp2);
}
