/* $Modified: Sat Jun  6 14:31:08 2009 by uwer $ */
#ifndef STANDARDMODELPARAMETERS_H
#define STANDARDMODELPARAMETERS_H


#include <map>
#include <iostream>
#include <cstdlib>
#include <limits>

enum PerturbativeOrder {LO = 0, NLO = 1, NNLO=2 };

enum Particles {
  GLUON = 0,
  DOWN=1, UP=2, STRANGE=3, CHARM=4, BOTTOM=5, TOP=6,
  ADOWN=7, AUP=8, ASTRANGE=9, ACHARM=10, ABOTTOM=11, ATOP=12, 
  ELECTRON=13, MUON=14, TAU=15, 
  WBOSON=16, ZBOSON=17, HIGGS=18, PHOTON=19,NUMBER_OF_PARTICLES=20
};


class StandardModelParameters {
 private:
  /*
   * We define StandardModelParameters as a Singleton:
   * The constructurs are all private, there is only one instance
   * of StandardModelParameters which is accessed through the memeber
   * function instance().
   * Some constants are definded as static const.
   * Note that the earlier version where StandardModelParameter was
   * not a Singleton but all the fields where static had problems
   * due to the unspecified initialization order.
   */

  double alphasLO;
  double alphasNLO;
  double FactorizationScale;
  double RenormalizationScale;
  double InternalScale;  
  double FactorizationScale2;
  double RenormalizationScale2;
  double InternalScale2;
  
  double Lambda5LO;
  double Lambda5NLO;
  double SinusThetaWeinbergSquared;
  double CosinusThetaWeinbergSquared;
  double alpha;
  
  double DeltaUV1;
  double DeltaIR1;
  double DeltaIR2;
  int ActiveFlavours;

  double masses[NUMBER_OF_PARTICLES];
  double masses2[NUMBER_OF_PARTICLES];
  std::map<Particles,double> AxialCouplings;
  std::map<Particles,double> VectorCouplings;
  std::map<Particles,double> charges;
  std::map<Particles,double> weakIsospins;
  

  double calculateAxialCoupling(Particles particle);
  double calculateVectorCoupling(Particles particle);
  void updateAlphas();
  void updateWeakCouplings();
  double AlphasEvolution(double mu, PerturbativeOrder order) const;

  /*
   * Make all these guys private to avoid generation of a 
   * Standardmodelobject 
   */ 
  StandardModelParameters();
  StandardModelParameters(const char* file, long line);
  StandardModelParameters(const StandardModelParameters &);
  void operator=(StandardModelParameters &);

 public:
  static const double hcq;
  static const double NumberOfColours;
  static const double CA;
  static const double CF;
  static const double TR;

  static inline StandardModelParameters & instance() {
    static StandardModelParameters standardmodelparameters;
    return standardmodelparameters;
  }

  void setFactorizationScale(double mu, int level = 0);
  void setRenormalizationScale(double mu, int level = 0);
  void setInternalScale(double mu,int level = 0);  

  void setMass(Particles particle,double mass);

  void setSinusThetaWeinbergSquared(double swq, int level = 0);

  void setAlpha(double val);

  void setLambda5LO(double lambda);
  void setLambda5LO(double scale, double alphas);

  static inline bool isquark(Particles t){
    if ( ( GLUON < t ) && (t <= ATOP ) )
      return true;
    return false;
  }

  static inline bool isgluon(Particles t){
    if ( GLUON == t )
      return true;
    return false;
  }

  static inline Particles antiquark(Particles t){
    if ( t < ADOWN) 
      return( static_cast<Particles>(t + ADOWN - DOWN) );
    if ( t > TOP ) 
      return( static_cast<Particles>(t + DOWN - ADOWN) );
    std::cout << "We should never end up here..." << std::endl;
    exit(1);
    return(GLUON);
  }

  inline double getHcq() const {
    return( hcq ) ;
  }

  inline double getNumberOfColours() const {
    return( NumberOfColours ) ;
  }

  inline double getSinusThetaWeinbergSquared() const {
    return( SinusThetaWeinbergSquared ) ;
  }

  inline double getCosinusThetaWeinbergSquared() const {
    return( CosinusThetaWeinbergSquared ) ;
  }

  inline double getAlpha() const {
    return( alpha ) ;
  }

  inline double getAxialCoupling(Particles particle)  {
    return( AxialCouplings[particle] ) ;
  }

  inline double getVectorCoupling(Particles particle) {
    return( VectorCouplings[particle] ) ;
  }
  
  inline double getFactorizationScale() const {
    return(FactorizationScale);
  }
  
  inline double getRenormalizationScale() const {
    return(RenormalizationScale);
  }
  
  inline double getInternalScale() const {
    return(InternalScale);
  }

  inline double getQaux2() const {
    return(RenormalizationScale2);
  }
  
  inline double getFactorizationScaleSquared() const {
    return(FactorizationScale2);
  }
  
  inline double getRenormalizationScaleSquared() const {
    return(RenormalizationScale2);
  }
  
  inline double getInternalScaleSquared() const {
    return(InternalScale2);
  }
 
  inline double getDeltaUV1() const {
    return(DeltaUV1);
  }
  
  inline double getDeltaIR1() const {
    return(DeltaIR1);
  }
  
  inline double getDeltaIR2() const {
    return(DeltaIR2);
  }
  
  inline void setDeltaUV1(double delta) {
    DeltaUV1=delta;
  }
  
  inline void setDeltaIR1(double delta) {
    DeltaIR1=delta;
  }
  
  inline void setDeltaIR2(double delta) {
    DeltaIR2=delta;
  }

  inline double getMass(Particles particle)  {  
    return(masses[particle]);
  }
  
  inline double getMassSquared(Particles particle)  {
    return(masses2[particle]);
  }

  // SA: Put explicitly a (quiet) NaN  to check that virtuals do
  // not depend on alphas
  inline double getAlphasLO() const{
    //return( alphasLO );
    return std::numeric_limits<double>::quiet_NaN();    
  }
  // SA: Put explicitly a (quiet) NaN  to check that virtuals do
  // not depend on alphas
  inline double getAlphasNLO() const{
    //return( alphasNLO );
    return std::numeric_limits<double>::quiet_NaN();    
  }

  void printParameters();
};


std::ostream& operator<<(std::ostream& os, const Particles& p);

#endif // STANDARDMODELPARAMETERS_H
