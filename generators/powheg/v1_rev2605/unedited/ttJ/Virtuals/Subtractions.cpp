// $Modified: Thu Jan 18 11:43:19 2007 by puwer $
#include "Subtractions.h"
#include "StandardModelParameters.h"
#include <cmath>
#include <vector>


namespace Subtractions 
{

  using namespace std;


  /***
   *** Interface to the Li2 fortran routine. Note that for x>1 only the
   *** real part is calculated.
   ***/
  extern "C" {
    double li2_(const double & x);
  }


  static StandardModelParameters&  parms = StandardModelParameters::instance();

  
  double Vsing(const double& sjk, const double& mj, const double& mk) {
    /*
     * We absorb the prefactor 
     *         1/Gamma(1-e)*(4pimu^2/sja)^e 
     * in the definition of Vsing. 
     * This allows us to use Delta1 and Delta2
     */ 
    
    // [hep-ph/0201036, eq. 6.20]
    if ( ( mj == 0.0 ) && ( mk == 0.0 ) ) {
      return( + parms.getDeltaIR2() 
	      + parms.getDeltaIR1() 
	      * log( parms.getInternalScaleSquared()/sjk)
	      + (1.0/2.0) * pow( log( parms.getInternalScaleSquared()/sjk),2) 
	      - pi2/6.0 );
    }
    
    if ( ( mj > 0.0 ) && ( mk == 0.0 ) ) {
      double mjq = mj * mj;
      double Qjk2 = sjk + mjq;
      return( + (1.0/2.0) * Vsing(sjk,0.0,0.0) 
	      + (1.0/2.0) * ( parms.getDeltaIR1() + 
			    log(parms.getInternalScaleSquared()/sjk) 
			    ) * log(mjq/sjk)
	      - (1.0/4.0) * pow(log(mjq/sjk),2) - pi2/12.0
	      - (1.0/2.0) * log(mjq/sjk)*log(sjk/Qjk2) 
	      - (1.0/2.0) * log(mjq/Qjk2)*log(sjk/Qjk2)
	      );
    }
    
    if ( ( mk > 0.0 ) && ( mj == 0.0 ) ) {
      double mkq = mk * mk;
      double Qjk2 = sjk + mkq;
      return( (1.0/2.0) * Vsing(sjk,0.0,0.0) 
	      + (1.0/2.0) * ( + parms.getDeltaIR1() 
			    + log(parms.getInternalScaleSquared()/sjk) 
			    ) * log(mkq/sjk)
	      - (1.0/4.0) * pow(log(mkq/sjk),2) - pi2/12.0
	      - (1.0/2.0) * log(mkq/sjk)*log(sjk/Qjk2) 
	      - (1.0/2.0) * log(mkq/Qjk2)*log(sjk/Qjk2)
	      );
    }

    if ( ( mj > 0.0 ) && ( mk > 0.0 ) ) {
      double mjq = mj * mj;
      double mkq = mk * mk;
      double Qjk2 = sjk + mjq + mkq;
      double vjk = sqrt( 1.0 - 4.0*mjq*mkq/sjk/sjk );
      double r = rho(vjk);
      double rj2 = rhoj2(vjk,mjq/Qjk2,mkq/Qjk2);
      double rk2 = rhok2(vjk,mjq/Qjk2,mkq/Qjk2);
      
      return( 1.0/vjk * ( + ( + parms.getDeltaIR1() 
			      + log(parms.getInternalScaleSquared()/sjk) 
			      ) * log(r)
			  - (1.0/4.0) * pow(log(rj2),2)
			  - (1.0/4.0) * pow(log(rk2),2) - pi2/6.0
			  + log(r) * log(Qjk2/sjk) ) );
    }
    
    cout << "Vsing: We should never end up here" << endl;
    cout << "Vsing: sjk = " << sjk << endl;
    cout << "Vsing: mk  = " << mk << endl;
    cout << "Vsing: mj  = " << mj << endl;
    exit(1);
  } 
  

  double VNSq(const double& sjk, const double& mj, const double& mk ) {
    
    if ( ( mj == 0.0 ) && (mk == 0.0 ) ) {
      return(0.0);
    }
    
    if ( ( mj > 0.0 ) && (mk == 0.0 ) ) {
      double mjq = mj * mj;
      double Qjk2 = sjk + mjq;
      // [hep-ph/0201036, eq. 6.22]
      return( gammaq/CF*log(sjk/Qjk2) + pi2/6.0 - li2_(sjk/Qjk2)
	      - 2.0*log(sjk/Qjk2) - mjq/sjk * log(mjq/Qjk2) );
    }
    
    if ( ( mj == 0.0 ) && (mk > 0.0 ) ) {
      double mkq = mk * mk;
      double Qjk2 = sjk + mkq;
      double Qjk = sqrt(Qjk2);
      // [hep-ph/0201036, eq. 6.23]
      return(gammaq/CF*(log(sjk/Qjk2) - 2.0*log( (Qjk-mk)/Qjk ) 
			- 2.0 * mk / ( Qjk + mk ) )
	     + pi2/6.0 - li2_(sjk/Qjk2) );
    }
    
    if ( ( mj > 0.0 ) && (mk > 0.0 ) ) {
      double mjq = mj * mj;
      double mkq = mk * mk;
      double Qjk2 = sjk + mjq + mkq;
      double Qjk = sqrt(Qjk2);
      double vjk = sqrt( 1.0 - 4.0*mjq*mkq/sjk/sjk );
      double r = rho(vjk);
      double r2 = r*r;
      double rj2 = rhoj2(vjk,mjq/Qjk2,mkq/Qjk2);
      double rk2 = rhok2(vjk,mjq/Qjk2,mkq/Qjk2);
      // [hep-ph/0201036, eq. 6.21]
      return(gammaq/CF*log(sjk/Qjk2) 
	     + 1.0/vjk * (log(r2)*log(1.0+r2) + 2.0 * li2_(r2) - li2_(1.0-rj2)
			  - li2_(1.0-rk2) - pi2/6.0 )
	     + log( (Qjk-mk)/Qjk ) - 2.0 * log( (pow(Qjk-mk,2)-mjq)/Qjk2 )
	     - 2.0*mjq/sjk * log( mj / (Qjk-mk) )
	     - mk / (Qjk-mk) + 2.0 * mk * ( 2.0*mk - Qjk ) / sjk + pi2/2.0 );
    }
    
    cout << "VNSq: We should never end up here" << endl;
    cout << "VNSq: sjk = " << sjk << endl;
    cout << "VNSq: mk  = " << mk << endl;
    cout << "VNSq: mj  = " << mj << endl;
    exit(1);
    
  }
  
  double VNSg(const double& sjk, const double& mj, const double& mk, 
	      vector<double>& mf, const double& kappa) {
    
    int nf = mf.size();
  
    vector<double> mfh;
   
    for (int i=0; i < nf; i++ ) {
      if ( sjk > 4.0 * mf[i]*(mf[i]+mk) ) {      
	mfh.push_back(mf[i]);
      }
    }
   
    if ( ( mj == 0.0 ) && ( mk == 0.0 ) ) {
      // [hep-ph/0201036, eq. 6.26]    
      if ( nf == 0 )
	return(0.0);
     
      double mkq = mk * mk;
      double Qjk2 = sjk + mkq;
      double Qjk = sqrt(Qjk2);
     
      double res = 0;
     
      for (vector<double>::iterator i = mfh.begin(); i!= mfh.end(); i++ ) {
	cout << "loop over massive Quarks...\n";
	double mfi = (*i);
	double r1 = sqrt(1.0 - 4.0 * mfi * mfi / pow(Qjk-mk,2) );
	double r2 = sqrt(1.0 - 4.0 * mfi * mfi/ ( Qjk2-mkq ) );
	res += (4.0/3.0)*TR/CA*( log( (1.0 + r1 )/ 2.0 ) - r1/3.0*(3.0+r1*r1)
			       - (1.0/2.0)*log(mfi*mfi/sjk) );
      }
      for (vector<double>::iterator i = mf.begin(); i!= mf.end(); i++ ) {
	cout << "loop over massive Quarks...\n";
	double mfi = (*i);
	res += (2.0/3.0)*TR/CA * log(mfi*mfi/parms.getQaux2());
      }
      return(res);
    }
   
    if ( ( mj == 0.0 ) && ( mk > 0 ) ) {
      double mkq = mk * mk;
      double Qjk2 = sjk + mkq;
      double Qjk = sqrt(Qjk2);
     
      double res = 
	+ gammag/CA * ( log(sjk/Qjk2) - 2.0 * log( (Qjk-mk) / Qjk )
			- 2.0 * mk / (Qjk+mk) ) 
	+ pi2/6.0 
	- li2_(sjk/Qjk2);
     
      for ( vector<double>::iterator i = mfh.begin(); i!= mfh.end(); i++ ) {
	cout << "loop over massive Quarks...\n";
	double mfi = (*i);
	double r1 = sqrt(1.0 - 4 * mfi / pow(Qjk-mk,2) );
	double r2 = sqrt(1.0 - 4 * mfi / ( Qjk2-mkq ) );
	res += (4.0/3.0) * TR / CA 
	  * ( log( (Qjk-mk) / Qjk ) + mk*pow(r1,3)/(Qjk+mk)
	      + log( (1.0+r1)/2.0 ) - r1/3.0*(3.0+r1*r1)
	      - 1.0/2.0*log(mfi*mfi/Qjk2) )
	  + ( kappa - (2.0/3.0) ) * mkq/sjk
	  * 2.0 * TR / CA * ( pow(r2,3)*log( (r2-r1)/(r2+r1) )
			      - log( ( 1.0-r1 )/( 1.0 + r1 ) )
			      - 8.0*r1*mfi*mfi/sjk );
      }
     
      for ( vector<double>::iterator i = mf.begin(); i!= mf.end(); i++ ) {
	cout << "loop over massive Quarks...\n";
	double mfi = (*i);
	res += (2.0/3.0)*TR/CA * log( mfi*mfi/parms.getQaux2() );
      }

      res += (kappa-2.0/3.0) * mkq/sjk * 
	2.0 * (TR/CA*NF-1.0)*log( 2.0 * mk / ( Qjk + mk ) );
     
      return(res);
     
    }
   
    cout << "VNSg: We should never end up here" << endl;
    cout << "VNSg: sjk = " << sjk << endl;
    cout << "VNSg: mk  = " << mk << endl;
    cout << "VNSg: mj  = " << mj << endl;
    exit(1);
   
  }


  /*
   * We absorb a factor 
   *
   * (4pi)^e/Gamma(1-e) 
   *
   * in the Gamma-Functions 
   */
  double Gamma_g(vector<double>& mf){
    double res = ( 
		  + parms.getDeltaIR1() 
		  + log( parms.getInternalScaleSquared()
			 / parms.getRenormalizationScaleSquared() ) 
		  ) * gammag;
    
    for ( vector<double>::iterator i = mf.begin(); i!= mf.end(); i++ ) {
      cout << "loop over massive Quarks...\n";
      double mfi = (*i);
      res -= (2.0/3.0)*TR * log( mfi*mfi/parms.getQaux2() );
    }
    return(res);
  }
  

  
}
