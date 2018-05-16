#include "StandardModelParameters.h"
#include <map>
#include <iostream>
#include <cmath>

using namespace std;


// 24.10.06 hcq=0.389 37966e-3 --> 0.389 379 323e-3
const double StandardModelParameters::hcq = 0.389379323e-3;
const double StandardModelParameters::NumberOfColours = 3.0;
const double StandardModelParameters::CA = 3.0;;
const double StandardModelParameters::CF = 4.0/3.0;
const double StandardModelParameters::TR = 0.5;


StandardModelParameters::StandardModelParameters(const char* file, long line){

#ifdef VERBOSE
  cout << "# -------------------------------------------------------------" 
       << endl;
  cout << "# Create instance of StandardModelParameters \n";
  cout << "# Called from file: " << file << " line: " << line << endl;
  cout << "# -------------------------------------------------------------" 
       << endl;
  cout << "#" << endl;
#endif

  StandardModelParameters();

}

StandardModelParameters::StandardModelParameters(){

  const double OneHalf  = (1.0/2.0);
  const double OneThird = (1.0/3.0);
  const double TwoThird = (2.0/3.0);

#ifdef VERBOSE
  cout << "#\n";
  cout << "# Initialize standard model parameters ..." << endl;
  
  cout << "# hcq = " << getHcq() << endl;
#endif
  Lambda5LO = 0.165;
  Lambda5NLO = 0.2260;
  
  alpha = 1.0/126.3;
  ActiveFlavours = 5;

  DeltaUV1 = 0.0;
  DeltaIR1 = 0.0;
  DeltaIR2 = 0.0;

  setMass( TOP, 174.0 );
  setMass( BOTTOM, 0. );
  setMass( CHARM, 0. );
  setMass( STRANGE, 0. );
  setMass( UP, 0.0 );
  setMass( DOWN, 0.0 );

  setMass( ELECTRON, 0. );
  setMass( MUON, 0. );
  setMass( TAU, 0. );
  
  setMass( WBOSON, 80.425 );
  setMass( ZBOSON, 91.1876 );
  setMass( HIGGS, 120.0 );
  
  setRenormalizationScale( masses[TOP]);
  setFactorizationScale( RenormalizationScale );
  setInternalScale( masses[TOP] );    
    
  /*
   * On-shell scheme for SinusThetaWeinberg
   */
  setSinusThetaWeinbergSquared(1.0 - masses2[WBOSON]/masses2[ZBOSON] );
  
  charges[TOP]   = TwoThird;
  charges[CHARM] = TwoThird;
  charges[UP]    = TwoThird;
  charges[BOTTOM]  = - OneThird;
  charges[STRANGE] = - OneThird;
  charges[DOWN]    = - OneThird;
  
  weakIsospins[TOP]   = OneHalf;
  weakIsospins[CHARM] = OneHalf;
  weakIsospins[UP]    = OneHalf;
  
  weakIsospins[BOTTOM]  = - OneHalf;
  weakIsospins[STRANGE] = - OneHalf;
  weakIsospins[DOWN]    = - OneHalf;
  
  updateWeakCouplings();
#ifdef VERBOSE
  cout << "# Standard model parameters initialized..." << endl;
#endif
  
}

void StandardModelParameters::setAlpha(double val){
  alpha = val;
  updateWeakCouplings();
}

void StandardModelParameters::updateAlphas(){
  alphasLO = AlphasEvolution(RenormalizationScale,LO);
  alphasNLO = AlphasEvolution(RenormalizationScale,NLO);
}

void StandardModelParameters::updateWeakCouplings(){

  AxialCouplings[TOP]     =  calculateAxialCoupling(TOP);
  AxialCouplings[BOTTOM]  =  calculateAxialCoupling(BOTTOM);
  AxialCouplings[CHARM]   =  calculateAxialCoupling(CHARM);
  AxialCouplings[STRANGE] =  calculateAxialCoupling(STRANGE);
  AxialCouplings[DOWN]    =  calculateAxialCoupling(DOWN);
  AxialCouplings[UP]      =  calculateAxialCoupling(UP);
  
  VectorCouplings[TOP]     =  calculateVectorCoupling(TOP);
  VectorCouplings[BOTTOM]  =  calculateVectorCoupling(BOTTOM);
  VectorCouplings[CHARM]   =  calculateVectorCoupling(CHARM);
  VectorCouplings[STRANGE] =  calculateVectorCoupling(STRANGE);
  VectorCouplings[DOWN]    =  calculateVectorCoupling(DOWN);
  VectorCouplings[UP]      =  calculateVectorCoupling(UP);
  
}

double StandardModelParameters::calculateAxialCoupling(Particles particle){
  return( 1.0 / 2.0 / sqrt(SinusThetaWeinbergSquared )
	  / sqrt(CosinusThetaWeinbergSquared ) * weakIsospins[particle] );
}

double StandardModelParameters::calculateVectorCoupling(Particles particle){
  return( 1.0 / 2.0 / sqrt(SinusThetaWeinbergSquared )
	  / sqrt(CosinusThetaWeinbergSquared ) 
	  * ( weakIsospins[particle] 
	      - 2.0 * SinusThetaWeinbergSquared * charges[particle] ) );
}

void StandardModelParameters::setFactorizationScale(double mu, int level){
  FactorizationScale = mu;
  FactorizationScale2 = FactorizationScale * FactorizationScale;
  if (level==1)
    cout << "#" << endl
	 << "# StandardModelParameters: Factorization scale = " << mu 
	 << endl << "#" << endl;
}

void StandardModelParameters::setRenormalizationScale(double mu, int level){
  RenormalizationScale = mu;
  RenormalizationScale2 = RenormalizationScale * RenormalizationScale;
  if (level==1)
    cout << "#" << endl
	 << "# StandardModelParameters: Renormalization scale = " << mu 
	 << endl << "#" << endl;
  updateAlphas();
}

void StandardModelParameters::setInternalScale(double mu, int level){
  InternalScale = mu;
  InternalScale2 = InternalScale * InternalScale;
  if (level==1)
    cout << "#" << endl
	 << "# StandardModelParameters: Internal scale = " << mu 
	 << endl << "#" << endl;
}

void StandardModelParameters::setSinusThetaWeinbergSquared(double swq, int level){
  SinusThetaWeinbergSquared = swq;
  CosinusThetaWeinbergSquared = 1.0 - swq;
  updateWeakCouplings();
  if (level==1) 
    cout << "#" << endl
	 << "# StandardModelParameters: SinusThetaWeinbergSquared = " 
	 << swq << endl
	 << "#" << endl;

}


void StandardModelParameters::setMass(Particles particle, double mass) 
{
  masses[particle] = mass;
  masses2[particle] = mass*mass;
}

void StandardModelParameters::setLambda5LO(double lambda){
  Lambda5LO = lambda;
  updateAlphas();
}

void StandardModelParameters::setLambda5LO(double scale, double alphas){
  const int nf = 5;
  const double beta0 = 11.0 - (2.0/3.0) * nf;

  Lambda5LO = scale / exp( ( 4. * M_PI ) / ( 2. * beta0 * alphas )  );

  cout << "# Set Lambda5LO =" << Lambda5LO << endl;

  updateAlphas();

  cout << "# alphasLO(" << scale << ") = " << getAlphasLO() << endl; 
}

double StandardModelParameters::AlphasEvolution(double mu, 
						PerturbativeOrder order) const
{ 
  const int nf = 5;
  const double beta0 = 11.0 - (2.0/3.0) * nf;
  const double beta1 = 51.0 - (19.0/3.0) * nf;
  const double beta2 = 2857.0 - (5033.0/9.0) * nf + (325.0/27.0)*nf*nf;

  if (order == LO) {
    return( 4.0 * M_PI / ( 2.0 * beta0 * log(mu/Lambda5LO) ) );
  }

  if (order == NLO) {
    const double lambda = Lambda5NLO;
    double alphas_lo = 4.0 * M_PI/( 2.0 * beta0 * log(mu/lambda) );
    double alphas_nlo = alphas_lo * 
      ( 1.0 - 2.0*beta1/beta0/beta0 * log(2.0 * log(mu/lambda) )
	/2.0/log(mu/lambda) ); 
    return(alphas_nlo);
  }
    
  cout << "Only LO and NLO implemented " 
       << __FILE__ << " " << __LINE__ << endl;
  exit(1);
}

void StandardModelParameters::printParameters(){
  cout << "#" <<endl;
  cout << "# -------------------------------------------------------------\n";
  cout << "#  Standard Model Parameters: " << endl;
  cout << "# -------------------------------------------------------------\n";
  cout << "#  Masses:" << endl;
  cout << "#                   mt = " << getMass(TOP) << endl;
  cout << "#                   mb = " << getMass(BOTTOM) << endl;
  cout << "#                   mc = " << getMass(CHARM) << endl;
  cout << "#                   ms = " << getMass(STRANGE) << endl;
  cout << "#                   mu = " << getMass(UP) << endl;
  cout << "#                   md = " << getMass(DOWN) << endl;
  cout << "#                   mw = " << getMass(WBOSON) << endl;
  cout << "#                   mz = " << getMass(ZBOSON) << endl;
  cout << "#                   mh = " << getMass(HIGGS) << endl;
  cout << "#" <<endl;
  cout << "#" <<endl;
  cout << "#  Couplings:" << endl;
  cout << "#                   alphasLO = " << getAlphasLO() << endl;
  cout << "#                  alphasNLO = " << getAlphasNLO() << endl;
  cout << "#                      alpha = " << getAlpha() << endl;
  cout << "#                    1/alpha = " << 1./getAlpha() << endl;
  cout << "#  SinusThetaWeinbergSquared = " << getSinusThetaWeinbergSquared() 
       << endl;
  cout << "#" << endl;
  cout << "#" << endl;
  cout << "#  Scales:"<<endl;
  cout << "#" << endl;
  cout << "#   FactorizationScale = " << getFactorizationScale() << endl;
  cout << "# RenormalizationScale = " << getRenormalizationScale() << endl;
  cout << "#        InternalScale = " << getInternalScale() << endl;
  cout << "#                Qaux  = " << sqrt(getQaux2()) << endl;
  cout << "#" << endl;
  cout << "#             DeltaUV1 = " << getDeltaUV1() << endl;
  cout << "#             DeltaIR1 = " << getDeltaIR1() << endl;
  cout << "#             DeltaIR2 = " << getDeltaIR2() << endl;
  cout << "#" << endl;
  cout << "# -------------------------------------------------------------\n";
  cout << "#" << endl;
  cout << "#" << endl;
}

ostream& operator << (ostream& os, const Particles& p){
  switch (p) {
  case GLUON:
    os << "g"; 
    break;
  case UP:
    os << "u"; 
    break;
  case DOWN:
    os << "d"; 
    break;
  case CHARM:
    os << "c"; 
    break;
  case STRANGE:
    os << "s"; 
    break;
  case BOTTOM:
    os << "b"; 
    break;
  case TOP:
    os << "t"; 
    break;
  case AUP:
    os << "u~"; 
    break;
  case ADOWN:
    os << "d~"; 
    break;
  case ACHARM:
    os << "c~"; 
    break;
  case ASTRANGE:
    os << "s~"; 
    break;
  case ABOTTOM:
    os << "b~"; 
    break;
  case ATOP:
    os << "t~"; 
    break;
  default:
    os << "unknow particle";
  }
  return os;
}

#ifdef MAIN
int main(){
  cout.precision(15);
  StandardModelParameters& parameters=StandardModelParameters::instance();
  cout << parameters.getMass(ATOP) << endl;
  parameters.setRenormalizationScale(parameters.getMass(TOP)/10.0);
  cout << parameters.getAlphasLO() << endl;
  cout << parameters.getAlphasNLO() << endl;
  parameters.printParameters();
}
#endif
