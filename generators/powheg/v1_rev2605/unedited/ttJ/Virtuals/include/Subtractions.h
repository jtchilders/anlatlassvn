/* $Modified: Fri Feb 16 16:44:38 2007 by puwer $ */
#ifndef SUBTRACTIONS_H_
#define SUBTRACTIONS_H_

#include <iostream>
#include <vector>
#include <cmath> 
#include "StandardModelParameters.h"


namespace Subtractions {

  static StandardModelParameters & _parms_ 
    = StandardModelParameters::instance();

  enum ftype { 
    DELTA=1, 
    REGULAR=2, 
    PLUSX=3, 
    PLUS1=4
  };

  const double pi2 = M_PI * M_PI;


  const double CA = 3.0;
  const double CF = 1.0/2.0/CA*(CA*CA-1.0);
  const double TR = 1.0/2.0;
  
  const double NF = 5.0;
  const double gammaq = 3.0/2.0 * CF;
  const double gammag = 11.0/6.0 * CA - 2.0/3.0 * TR * NF;

  // [hep-ph/0201036, eq. 6.17]
  const double Kq = (7.0/2.0 - pi2/6.0) * CF;
  const double Kg = (67.0/18.0 - pi2/6.0) * CA - 10.0/9.0 * TR * NF;


  inline double Pprime_qq(const double& x){
    return( CF * (1.0-x) );
  }

  inline double Pprime_qg(const double& x){
    return( CF * x );
  }

  inline double Pprime_gq(const double& x){
    return( 2.0 * TR * x * (1.0 - x) );
  }

  inline double Pprime_gg(const double& x){
    return( 0.0 );
  }

  inline double Pqq_reg(const double& z) {
    return( -CF * ( 1.0 + z ) );
  }
  
  inline double Pgq_reg(const double& z) {
    return( TR * ( z*z + (1.0-z)*(1.0-z) ) );
  }

  inline double Pqg_reg(const double& z) {
    return( CF * (1.0 + (1.0-z) * (1.0-z) ) / z );
  }

  inline double Pgg_reg(const double& z) {
    return( 2.0 * CA * ( ( 1.0-z ) / z - 1.0 + z * (1.0 - z ) ) );
  }
  
  inline double Pgg_nonreg(const double& z) {
    return( 2.0 * CA / (1.0-z) );
  }
  
  inline double Pqq_nonreg(const double& z) {
    return( 2.0 * CF / (1.0-z) );
  }

  /*
   * The Altarelli-Paris functions, p.31 [hep-ph/0201036]
   */
  // [hep-ph/0201036, eq. 5.89]

  inline double Pgg(const double& x, const ftype& t) {

    if ( t == REGULAR ) 
      return( Pgg_reg(x) );

    if ( t == DELTA ) 
      return(gammag);

    if ( ( t == PLUSX ) || ( t == PLUS1 ) ) 
      return( Pgg_nonreg(x) );

    std::cout << "Pgg: We should never end up here...\n";
    exit(1);
  }

  inline double Pqq(const double& x, const ftype& t) {

    if ( t == REGULAR ) 
      return(Pqq_reg(x));

    if ( t == DELTA ) 
      return(gammaq);

    if ( ( t == PLUSX ) || ( t == PLUS1 ) ) 
      return( Pqq_nonreg(x) );

    std::cout << "Pqq: We should never end up here...\n";
    exit(1);

  }

  inline double Pqg(const double& x, const ftype& t) {

    if ( t == REGULAR ) 
      return(Pqg_reg(x));

    if ( ( t == DELTA ) || ( t == PLUSX ) || ( t == PLUS1 ) ) 
      return( 0.0 ); 

    std::cout << "Pqg: We should never end up here...\n";
    exit(1);

  }

  inline double Pgq(const double& x, const ftype& t) {

    if ( t == REGULAR ) 
      return(Pgq_reg(x));

    if ( ( t == DELTA ) || ( t == PLUSX ) || ( t == PLUS1 ) ) 
      return( 0.0 ); 

    std::cout << "Pgq: We should never end up here...\n";
    exit(1);

  }

  inline double rho(const double& vijk) {
    // [hep-ph/0201036, eq. 5.30]
    return(sqrt( ( 1.0 - vijk) / ( 1.0 + vijk) ) );
  }
  
  inline double rhoj2(const double& vijk, const double& mujq, 
		      const double& mukq) {
    // [hep-ph/0201036, eq. 5.30]
    return( ( 1.0 - vijk + 2.0 * mujq / (1.0 - mujq - mukq ) ) 
	    / ( 1.0 + vijk + 2.0 * mujq / (1.0  - mujq - mukq ) )  );
  }
  
  inline double rhok2(const double& vijk, const double& mujq, 
		      const double& mukq) {
    // [hep-ph/0201036, eq. 5.30]
    return( ( 1.0 - vijk + 2.0 * mukq / (1.0 - mujq - mukq ) ) 
	    / ( 1.0 + vijk + 2.0 * mukq / (1.0  - mujq - mukq ) ) );
  }

  double Vsing(const double& sjk, const double& mj, const double& mk);

  double VNSq(const double& sjk, const double& mj, const double& mk );

  double VNSg(const double& sjk, const double& mj, const double& mk, 
	      std::vector<double>& mf, const double& kappa);

  inline double Vq(const double& sjk, const double& mj, const double& mk ){
    return( Vsing(sjk,mj,mk) + VNSq(sjk,mj,mk) );
  }

  inline double Vg(const double& sjk, const double& mj, const double& mk, 
	    std::vector<double>& mf, const double& kappa){
    return( Vsing(sjk,mj,mk) + VNSg(sjk,mj,mk,mf,kappa) );
  }

  double Gamma_g(std::vector<double>& mf);



  inline double JgQ(const double& x, const double& muq, const ftype& t){

    /*
     * [hep-ph/0201036] eq.(5.58)
     */

    if ( ( t == DELTA ) || ( t == REGULAR ) )
      return 0.0;

    if ( t == PLUSX ) 
      return( ( + (1.0 - x) / 2.0 / pow(1.0 - x + muq,2) 
		- 2.0 / (1.0-x) * ( 1.0 + log( 1.0 - x + muq ) ) )
	      + 2.0 / ( 1.0 - x) * log( 2.0 + muq - x ) ); 

    if ( t == PLUS1 ) 
      return( ( + (1.0 - x) / 2.0 / pow(1.0 - x + muq,2) 
		- 2.0/(1.0-x) * ( 1.0 + log( 1.0 - x + muq ) ) )
	      + 2.0 / ( 1.0 - x) * log( 1.0 + muq ) ); 

    std::cout << "JgQ: We should never end up here...\n";
    exit(1);
  }

  inline double Kbargg(const double& x, const ftype& t){

    /*
     * [hep-ph/0201036] eq.(6.56)
     */

    if ( t == REGULAR ) {
      return( Pgg_reg(x) * log( ( 1.0 - x ) / x ) + Pprime_gg(x) ); 
    } 

    if ( t == DELTA ) {
      return( - ( gammag + Kg - 5.0/6.0 * pi2 * CA ) );
    }

    if ( ( t == PLUSX ) || ( t == PLUS1 )  ) {
      return( CA * 2.0 / ( 1.0 - x ) * log( ( 1.0 - x ) / x ) );
    }
    std::cout << "Kbargg: We should never end up here...\n";
    exit(1);
  }  

  inline double Kbarqq(const double& x, const ftype& t){

    /*
     * [hep-ph/0201036] eq.(6.56)
     */

    if ( t == REGULAR ) {
      return( Pqq_reg(x)*log( ( 1.0 - x ) / x ) + Pprime_qq(x) ); 
    }  
    if ( t == DELTA ) {
      return(- ( gammaq + Kq - 5.0/6.0 * pi2 * CF ) );
    }

    if ( ( t == PLUSX ) || ( t == PLUS1 )  ) {
      return( CF * 2.0 / ( 1.0 - x ) * log( ( 1.0-x ) /x ) );
    }
    std::cout << "Kbarqq: We should never end up here...\n";
    exit(1);
  }

  inline double Kbarqg(const double& x, const ftype& t){

    /*
     * [hep-ph/0201036] eq.(6.56)
     */

    if ( t == REGULAR ) {
      return( Pqg_reg(x) * log( ( 1.0 - x ) / x ) + Pprime_qg(x) ); 
    } 

    if ( ( t == DELTA ) || ( t == PLUSX ) || ( t == PLUS1 )  ) {
      return( 0.0 );
    }
    std::cout << "Kbarqg: We should never end up here...\n";
    exit(1);

  }

  inline double Kbargq(const double& x, const ftype& t){

    /*
     * [hep-ph/0201036] eq.(6.56)
     */

    if ( t == REGULAR ) {
      return( Pgq_reg(x)*log( (1.0-x)/x ) + Pprime_gq(x) ) ;
    } 

    if ( ( t == DELTA ) || ( t == PLUSX ) || ( t == PLUS1 )  ) {
      return( 0.0 );
    }
    std::cout << "Kbargq: We should never end up here...\n";
    exit(1);
  }

  inline double KRgg(const double& x, const double& sja, const double& mjq, 
	      const ftype& t){
    /*
     * This is the square bracket in [hep-ph/0201036] eq. 6.55
     */

    if ( mjq == 0. )
      return(0.);

    if ( t == REGULAR ){
      return( Pgg_reg(x) * log( (1.0-x)*sja / ( (1.0-x)*sja + mjq ) ));
    }

    if ( t == DELTA ){
      double mj = sqrt(mjq); 
      return(gammag * ( log( (sja - 2.0 * mj * sqrt(sja+mjq) + 2.0 * mjq ) 
			     / sja  ) + 2.0 * mj / ( sqrt(sja+mjq) + mj ) ) );
    }

    if ( ( t==PLUSX ) || (t==PLUS1) ){
      return(0.0);
    }

    std::cout << "KRgg: We should never end up here...\n";
    exit(1);
  }

  inline double KRqq(const double& x, const double& sja, const double& mjq, 
	      const ftype& t){
    /*
     * This is the square bracket in [hep-ph/0201036] eq. 6.55
     */

    if ( mjq == 0. )
      return(0.);

    if ( t == REGULAR ){
      return( Pqq_reg(x) * log( (1.0-x)*sja / ( (1.0-x)*sja + mjq ) ));
    }

    if ( t == DELTA ){
      double mj = sqrt(mjq); 
      return(gammaq * ( log( (sja - 2.0 * mj * sqrt(sja+mjq) + 2.0 * mjq ) 
			     / sja  ) + 2.0 * mj / ( sqrt(sja+mjq) + mj ) ) );
    }

    if ( ( t==PLUSX ) || (t==PLUS1) ){
      return(0.0);
    }

    std::cout << "KRqq: We should never end up here...\n";
    exit(1);
  }

  inline double KRqg(const double& x, const double& sja, const double& mjq, 
	      const ftype& t){
    /*
     * This is the square bracket in [hep-ph/0201036] eq. 6.55
     */

    if ( mjq == 0. )
      return(0.);

    if ( t == REGULAR ){
      return( Pqg_reg(x) * log( (1.0-x)*sja / ( (1.0-x)*sja + mjq ) ));
    }


    if ( ( t == DELTA ) || ( t==PLUSX ) || (t==PLUS1) ){
      return(0.0);
    }

    std::cout << "KRqg: We should never end up here...\n";
    exit(1);
  }

  inline double KRgq(const double& x, const double& sja, const double& mjq, 
	      const ftype& t){
    /*
     * This is the square bracket in [hep-ph/0201036] eq. 6.55
     */

    if ( mjq == 0. )
      return(0.);

    if ( t == REGULAR ){
      return( Pgq_reg(x) * log( (1.0-x)*sja / ( (1.0-x)*sja + mjq ) ));
    }


    if ( ( t == DELTA ) || ( t==PLUSX ) || (t==PLUS1) ){
      return(0.0);
    }

    std::cout << "KRgq: We should never end up here...\n";
    exit(1);
  }

  inline double calK_q_gq(const double& , const double& , const double& , 
		   const ftype& ){
    /*
     * [hep-ph/0201036] eq.(6.57)
     */
    return(0.);
  }

  inline double calK_g_gq(const double&, const double&,
			  std::vector<double>&, const ftype&){
    /*
     * [hep-ph/0201036] eq.(6.61)
     */
    return(0.);
  }

  inline double calK_g_qg(const double& x, const double& sja,
			  std::vector<double>& mf, const ftype& t){
    /*
     * [hep-ph/0201036] eq.(6.61)
     */
    return(0.);
  }

  inline double calK_q_qq(const double& x, const double& sja, 
			  const double& mjq, 
			  const ftype& t){

    /*
     * [hep-ph/0201036] eq.(6.58)
     */
    if ( t == REGULAR ) {
      if ( mjq == 0.0 ) 
	return(0.0);
      else
	return( - 2.0 * log(2.0-x)/( 1.0 - x ) ) ;
    } 

    if ( t == DELTA ) {
      if ( mjq == 0.0 ) 
	return( -gammaq/CF );
      else 
	return( -gammaq/CF 
		+ mjq/sja*log( mjq/(sja+mjq) ) + 1.0/2.0*mjq/(sja+mjq) );
    }

    if ( t == PLUSX ) {
      if ( mjq == 0.0 ) 
	return( - gammaq/CF / (1.0-x) ) ;
      else
	return( + 2.0*log(1.0-x)/(1.0-x) + JgQ(x,mjq/sja,PLUSX) 
		+ 2.0/(1.0-x)*log( (2.0-x)*sja/ ( (2.0-x)*sja+mjq ) ) );
    }

    if ( t == PLUS1 ) {
      if ( mjq == 0.0 ) 
	return( - gammaq/CF / (1.0-x) ) ;
      else
	return( + 2.0*log(1.0-x)/(1.0-x) + JgQ(x,mjq/sja,PLUS1) 
		+ 2.0/(1.0-x)*log( sja / ( sja+mjq ) ) );
    }

    std::cout << "calK_q_qq: We should never end up here...\n";
    exit(1);
  }

  inline double calK_q_qg(const double& x, const double& sja, 
			  const double& mjq, 
			  const ftype& t){

    /*
     * [hep-ph/0201036] eq.(6.59)
     */

    if ( t == REGULAR ) {
      if (mjq == 0.0) 
	return(0.0);
      else
	return( 2.0 * CF / CA * mjq / x / sja 
		* log( mjq / ( (1.0-x)*sja + mjq ) ) ); 
    } 

    if ( t == DELTA ) {
      return( 0.0 );
    }

    if ( t == PLUSX ) {
      return( 0.0 );
    }

    if ( t == PLUS1 ) {
      return( 0.0 );
    }

    std::cout << "calK_q_qg: We should never end up here...\n";
    exit(1);
  }


  inline double calK_g_gg(const double& x, const double& sja, 
		   std::vector<double>& mf, 
		   const ftype& t){
    /*
     * [hep-ph/0201036] eq.(6.61)
     */
    if ( mf.size() > 0 ) {
      std::cout << "calK_g_gg: mf.size() > 0 not yet implemented\n";
      exit(1);
    }
    if ( t == REGULAR ) 
      return( 0.0 );

    if ( t == DELTA ) 
      return( -gammag/CA );

    if ( ( t == PLUSX ) || ( t == PLUS1 ) ) 
      return( - gammag / CA / ( 1.0 - x ) ); 

    std::cout << "calK_g_gg: We should never end up here...\n";
    exit(1);
  }

  inline double calK_g_qq(const double& x, const double& sja, 
		   std::vector<double>& mf, 
		   const ftype& t){
    /*
     * [hep-ph/0201036] eq.(6.61)
     */
    if ( mf.size() > 0 ) {
      std::cout << "calK_g_qq: mf.size() > 0 not yet implemented\n";
      exit(1);
    }
    if ( t == REGULAR ) 
      return( 0.0 );

    if ( t == DELTA ) 
      return( -gammag/CA );

    if ( ( t == PLUSX ) || ( t == PLUS1 ) ) 
      return( - gammag / CA / ( 1.0 - x ) ); 

    std::cout << "calK_g_qq: We should never end up here...\n";
    exit(1);
  }

  inline double calK_q_gg(const double& x, const double& sja, 
			  const double& mjq, const ftype& t){
    /*
     * [hep-ph/0201036] eq.(6.60)
     */

    return( calK_q_qq(x,sja,mjq,t) + CA/CF * calK_q_qg(x,sja,mjq,t) );
  }

  inline double Gamma_g(){
    return( ( + _parms_.getDeltaIR1() 
	      + log(_parms_.getInternalScaleSquared()
		    /_parms_.getRenormalizationScaleSquared()) ) * gammag );
  }
  
  inline double Gamma_q(){
    return(( + _parms_.getDeltaIR1() 
	     + log(_parms_.getInternalScaleSquared()
		   /_parms_.getRenormalizationScaleSquared()) ) * gammaq);
  }
  
  inline double Gamma_q(const double& m){
    return( CF * ( +  
		   ( + _parms_.getDeltaIR1() 
		     + log(_parms_.getInternalScaleSquared()
			   /_parms_.getRenormalizationScaleSquared()) 
		     ) 
		   + 1.0/2.0 * log(m*m
				   /_parms_.getRenormalizationScaleSquared()) 
		   - 2.0 
		   ) 
	    );
  }

}

#endif // SUBTRACTIONS_H_
