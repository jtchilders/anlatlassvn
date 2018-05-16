#include <iostream>
#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <sstream>
#include "FourMomentum.h"
#include "SquaredMatrixElements.h" 
#include "FiniteVirtuals.h"
#include "GGReduction.h"
#include "SysInfo.h"
#include "Errors.h"
#include "EvalSMExx.h"
#include "Subtractions.h"
#include "Color.h"
#include "Reduction.h"
#include "ScalarInt.h"
#include "Evalf.h"
#include "StandardModelParameters.h"
#include <limits>
#define PR(xxx_) cout << #xxx_ << " = " << xxx_<< endl;
using namespace std;


static SysInfo sysinfo;

static StandardModelParameters& parms = StandardModelParameters::instance();

bool CheckLargeCorrections=true;
bool FFTestFlag=false;
double LargeCorrFact=100.0;


void printBanner(){

  cout<<"#"<<endl;
  cout<<"#"<<endl;
  cout<<"# Initialization of library for virtual routines."<<endl; 
  cout<<"# If you use this program towards a scientific publication,"<<endl; 
  cout<<"# please opportunely cite the original papers where the routines "<<endl;
  cout<<"# included in this library have been presented: "<<endl;
  cout<<"#"<<endl;
  cout<<"# S. Dittmaier, P. Uwer and S. Weinzierl,"<<endl;
  cout<<"# NLO QCD corrections to t anti-t + jet production at hadron colliders,"<<endl;
  cout<<"# Phys. Rev. Lett. 98, 262002 (2007)"<<endl;
  cout<<"# [hep-ph/0703120 [HEP-PH]]."<<endl;
  cout<<"#"<<endl;
  cout<<"# S. Dittmaier, P. Uwer and S. Weinzierl,"<<endl;
  cout<<"# Hadronic top-quark pair production in association with a hard jet at next to"<<endl;
  cout<<"# leading order QCD: Phenomenological studies for the Tevatron and the LHC,"<<endl;
  cout<<"# Eur. Phys. J. C  59, 625 (2009)"<<endl;
  cout<<"# [arXiv:0810.0452 [hep-ph]]."<<endl;
  cout<<"#"<<endl;
  cout<<"#"<<endl;

}

void initVirtualRoutines(){

  // Assign a given kinematic and call all the routines once at the
  // beginning, in order to perform the initialization (and printout)
  // steps here, before entering POWHEG routines

  FourMomentum p_1(458.533175385278300,207.025516990944000,
		   370.293273289616700,0.000000000000000);
  
  FourMomentum p_2(-500.000000000000000,0.000000000000000,
		   500.000000000000000,0.000000000000000);
  
  FourMomentum p_3(-500.000000000000000,0.000000000000000,
		   -500.000000000000000,0.000000000000000);

  FourMomentum p_4(334.866822006721700,-196.368580218418100,
		    -267.893452247508300,-42.523727809261470);
  
  FourMomentum p_5(206.600002608000000,-10.656936772525890,
		    -102.399821042108500,42.523727809261470);
  

  double tmp;
  cout<<"#"<<endl;

  tmp=gggtt2virtfin(-p_2,-p_3,p_4,p_1,p_5); 
  tmp=qgttq2virtfin(-p_2,-p_3,p_4,p_1,p_5); 
  tmp=gqbttqb2virtfin(-p_2,-p_3,p_4,p_1,p_5); 
  tmp=qqttg2virtfin(-p_2,-p_3,p_4,p_1,p_5); 

  cout<<"#"<<endl;
  cout<<"# Initialization of virtual library completed ! "<<endl; 
  cout<<"#"<<endl;

}


void SetVirtualScales(int verbose = 0);
void VirtualInitialize();

/// F77 to C++ wrappers

extern "C" {   

  void callback_(){};
  
  int f77spaversion_();

    
  // pwhg_st common block in POWHEG-BOX/include/pwhg_st.h
  extern struct {
    double st_alpha,st_mufact2,st_muren2,st_lambda5MSB, st_facfact,st_renfact;
    int st_nlight ;
  } pwhg_st_;


    // ph_common common block in PhysPars.h
  extern struct {
    double ph_alphaem,ph_Zmass,ph_Zwidth,ph_Wmass,ph_Wwidth,ph_cthw,ph_sthw,ph_sthw2,ph_topmass,ph_topmass2,ph_topwidth,ph_topmass2low,ph_topmass2high, ph_topmtopw,ph_unit_e,ph_CKM;
  } ph_common_;


  // ccheckvirtuals common block in POWHEG-BOX/ttJ/virtual.f
  extern struct {
    double bcut;
    double largecorrfact;
    bool checklargecorr;
  } ccheckvirtuals_;

  // cfftestflag common block in POWHEG-BOX/ttJ/virtual.f
  // Kept separate from ccheckvirtuals because LOGICAL in fortran // are
  // KIND=4, while bool in C/C++ are just one bit. The results is that
  // the alignment of common block is ruined after the inclusion of a
  // logical variable
  extern struct {
    bool fftestflag;
  } cfftestflag_;



  void virtual_initialize_(){
    VirtualInitialize();
  }   

 
  void set_virtual_scales_(){
    // Don't assign logic from fortran to bool in C++ : SIZE MATTERS 
    static bool firstcall=true;
    if((!ccheckvirtuals_.checklargecorr) && firstcall) {
      CheckLargeCorrections=false;
      cout<<RED<<"Large Correction check switched off for testing"<<RESET<<endl;
    }
    SetVirtualScales();
    firstcall=false;
  }

 
 // f77 interfaces to StandarModelParameters.h routines 
  void setdeltair2_(const double & value){ 
    parms.setDeltaIR2(value);
    //   cout<<RED<< "#DeltaIR2 = "<<parms.getDeltaIR2()<<RESET<<endl;
  }
  void setdeltair1_(const double & value){ 
    parms.setDeltaIR1(value);
    // cout<<RED<< "#DeltaIR1 = "<<parms.getDeltaIR1()<<RESET<<endl;
  }
  void setdeltauv1_(const double & value){ 
    parms.setDeltaUV1(value);
    // cout<<RED<< "#DeltaUV1 = "<<parms.getDeltaUV1()<<RESET<<endl;
  }
}

void VirtualInitialize() {

  printBanner();
  
  // Initialization of virtuals routines
  cout<<"#"<<endl;
  cout<<"# Runtime assignment of parameters in virtual library " <<endl;
  cout<<"#"<<endl;

  // Set the values of SM parameters in StandarModelParameters
  // parms is a singleton
  parms.setMass(TOP,ph_common_.ph_topmass);
  // Following quantities are irrelevant for this process!
  // but set them equal to POWHEG ones anyway
  parms.setMass( WBOSON, ph_common_.ph_Wmass);
  parms.setMass( ZBOSON, ph_common_.ph_Zmass );
  parms.setSinusThetaWeinbergSquared(ph_common_.ph_sthw2);
  parms.setAlpha(ph_common_.ph_alphaem);
  parms.setLambda5LO(pwhg_st_.st_lambda5MSB);
  
  SetVirtualScales();

  // Printout parameters

  parms.printParameters();

  //Set cutoff in Reduction

  if (ccheckvirtuals_.bcut<0.0) {
    setBCUT(1.e-14);
  } else {
    setBCUT(ccheckvirtuals_.bcut);
  }
      
  if (ccheckvirtuals_.largecorrfact<0.0) {
    cout << "# Set LargeCorrFact = " << LargeCorrFact << endl;
  } else {
    LargeCorrFact=ccheckvirtuals_.largecorrfact;
    cout << "# Set LargeCorrFact = " << LargeCorrFact << endl;
  }


  //Some checks
    
  if(parms.getDeltaIR2() != 0) {
    cout<<RED<< "# ERROR DeltaIR2 = "<<parms.getDeltaIR2()<<RESET<<endl;
    exit(1);
  } 

  if(parms.getDeltaIR1() != 0) {
    cout<<RED<< "# ERROR DeltaIR1 = "<<parms.getDeltaIR1()<<RESET<<endl;
    exit(1);
  } 
    
  if(parms.getDeltaUV1() != 0) {
    cout<<RED<< "# ERROR DeltaUV1 = "<<parms.getDeltaUV1()<<RESET<<endl;
    exit(1);
  } 


  //Check the definitions of spinor products

  if (0==spaversion()){
    cout << "#\n";
    cout << "# Using the old definition for the spinor products\n";
    cout << "# singular for momenta along the z-axis\n";
  } else {
    cout << "#\n";
    cout << "# Using the new definition for the spinor products\n";
  }
  
  if (spaversion() != f77spaversion_()){
    cout << "# Inconsistent built: spinor product version used in f77 part\n";
    cout << "# does not match C++ version, rebuilt program!\n";
    exit(1);
  } else {
    cout << "# Checking consistency of spinor product version used "
	 << "in f77 part...\n";
    cout << "# Same version!"<<endl;
  }

  // Initialize  integral library
 
   if(cfftestflag_.fftestflag) {
      FFTestFlag=true;
  }

  initIntegralLibrary();

  // Call the virtual routines in order to initialize the cache

  initVirtualRoutines();

}


void SetVirtualScales(int verbose) {
  // Assign scales and check alpha_s on a event-by-event basis
  
  double mur = sqrt(pwhg_st_.st_muren2); 

  // Put explicitly a (quiet) NaN in muf to check that virtuals do
  // not depend on the factorization scale 
  // double muf = sqrt(pwhg_st_.st_mufact2);
  double muf=numeric_limits<double>::quiet_NaN();    
    
  double mt  = ph_common_.ph_topmass;

  // the second argument below is the verbosity level
  // 0: no output (default)
  // 1: verbose 
  parms.setInternalScale(mur,verbose);
  parms.setRenormalizationScale(mur,verbose);
  parms.setFactorizationScale(muf,verbose);
   
  /*SA: Commented out to allow for different PDF sets, with their own alphas definition
    if (abs(pwhg_st_.st_alpha - parms.getAlphasNLO())/pwhg_st_.st_alpha > 1.e-14) {
    cout.precision(16);
    cout<<RED<< "ERROR in alphas value assignment in virtual_initialize_() "<<RESET<<endl;
    cout<<RED<< pwhg_st_.st_alpha << " "<< parms.getAlphasNLO() <<RESET<<endl;
    cout<<RED<< pwhg_st_.st_muren2 << " "<< mur*mur <<RESET<<endl;
    exit(1);
    };*/
}


extern "C" {   

  // f77 interface to modified virtuals written by SA
  //
  // Corresponds to the following Fortran subroutine
  // gggtt2virtfin(p(0,1),p(0,2),p(0,5),p(0,3),p(0,4))
  //
  double gggtt2virtfin_(const double *p1,const double *p2, const double *p3,const double *kq,const double *kqb){
    FourMomentum xp1(p1), xp2(p2), xp3(p3), xkq(kq), xkqb(kqb);

    double res,res2;
    
    res=gggtt2virtfin(xp1,xp2,xp3,xkq,xkqb);
   
//     res2=gggtt2virtfin(xp1,xp2,xp3,xkq,xkqb);
//     if ((abs(res-res2)/res)>1e-7) {
//       cout.precision(16);
//       cout<<RED<<res<<" "<<res2<<RESET<<endl;
//        for (int i=0;i<10;i++) {
// 	 cout<<RED<<gggtt2virtfin(xp1,xp2,xp3,xkq,xkqb)<<RESET<<endl;
//        }
//        exit(1);
//     }
    return res;
  }

  // f77 interface to modified virtuals written by SA
  //
  // Corresponds to the following Fortran subroutine
  // gqbttqb2virtfin(p(0,1),p(0,2),p(0,5),p(0,3),p(0,4))
  //
  double gqbttqb2virtfin_(const double *p1,const double *p2, const double *p3,const double *kq,const double *kqb){
    FourMomentum xp1(p1), xp2(p2), xp3(p3), xkq(kq), xkqb(kqb);

    double res;
    

//#define BADCONFIGURATION
#ifdef BADCONFIGURATION    
    FourMomentum p_1(947.295708421587619,947.295708421587619,
		     0.000000000000000,0.000000000000000);
    
    FourMomentum p_2(947.295708421587619,-947.295708421587619,
		     0.000000000000000,0.000000000000000);
  
    FourMomentum p_3(943.316635236800039,-2.856123836693827,
		     -660.624088762582574,-650.490530459679235);

    FourMomentum p_4(946.217751101886165,-0.044992248263762,
		      664.766199696023818,650.490530459679235);
  
    FourMomentum p_5(5.057030504489035,2.901116084957589,
		     -4.142110933441247,0.000000000000000);



    cout<<"kt p3 p4  "<<p_3.kt(p_4)<<endl;
    cout<<"kt p4 p3  "<<p_4.kt(p_3)<<endl;
    cout<<"kt p4 p5  "<<p_4.kt(p_5)<<endl;
    cout<<"kt p5 p4  "<<p_5.kt(p_4)<<endl;
    cout<<"kt p3 p5  "<<p_3.kt(p_5)<<endl;
    cout<<"kt p5 p3  "<<p_5.kt(p_3)<<endl;
    
    cout<<"p3^2 "<<dotp(p_3,p_3)<<endl;
    cout<<"p4^2 "<<dotp(p_4,p_4)<<endl;
    cout<<"p5^2 "<<dotp(p_5,p_5)<<endl;
    cout<<"finmom "<<p_3+p_4+p_5<<endl;
    
    res=gqbttqb2virtfin(p_2,p_1,p_5,p_3,p_4);

    cout<<RED<<res<<RESET<<endl;

#else
    //xp1.printMaple();
    //xp2.printMaple();
    //xp3.printMaple();
    //xkq.printMaple();
    //xkqb.printMaple();
    res=gqbttqb2virtfin(xp1,xp2,xp3,xkq,xkqb);
    //exit(1);
#endif
  
    //cout<<RED<<res<<RESET<<endl;

    return res;
  }
 

  // f77 interface to modified virtuals written by SA
  //
  // Corresponds to the following Fortran subroutine
  // qgttq2virtfin(p(0,1),p(0,2),p(0,5),p(0,3),p(0,4))
  //
  double qgttq2virtfin_(const double *p1,const double *p2, const double *p3,const double *kq,const double *kqb){
    FourMomentum xp1(p1), xp2(p2), xp3(p3), xkq(kq), xkqb(kqb);

    double res;
    
    res=qgttq2virtfin(xp1,xp2,xp3,xkq,xkqb);
  
    //  cout<<RED<<res<<RESET<<endl;

    return res;
  }


  // f77 interface to modified virtuals written by SA
  //
  // Corresponds to the following Fortran subroutine
  // qqttg2virtfin(p(0,1),p(0,2),p(0,5),p(0,3),p(0,4))
  //
  double qqttg2virtfin_(const double *p1,const double *p2, const double *p3,const double *kq,const double *kqb){
    FourMomentum xp1(p1), xp2(p2), xp3(p3), xkq(kq), xkqb(kqb);

    double res;
    
    res=qqttg2virtfin(xp1,xp2,xp3,xkq,xkqb);
  
    //  cout<<RED<<res<<RESET<<endl;

    return res;
  }
 

  /////////////////////////////////////////////////////////////
  ///    THE FOLLOWING ROUUTINES ARE USED ONLY FOR CHECKS   ///
  /////////////////////////////////////////////////////////////

  // f77 interface to Born in Uwer's library
  //
  // Corresponds to the following Fortran subroutine
  // gggtt2(p(0,1),p(0,2),p(0,5),p(0,3),p(0,4))
  //
  double gggtt2_(const double *p1,const double *p2, const double *p3,const double *kq,const double *kqb){
    FourMomentum xp1(p1), xp2(p2), xp3(p3), xkq(kq), xkqb(kqb);
    double res;

    res=gggtt2(xp1,xp2,xp3,xkq,xkqb);
  
    // cout<<RED<<res<<RESET<<endl;

    return res;
  }


  // f77 interface to virtuals in Uwer's library
  //
  // Corresponds to the following Fortran subroutine
  // gggtt2virt(p(0,1),p(0,2),p(0,5),p(0,3),p(0,4))
  //
  double gggtt2virt_(const double *p1,const double *p2, const double *p3,const double *kq,const double *kqb){
    FourMomentum xp1(p1), xp2(p2), xp3(p3), xkq(kq), xkqb(kqb);

    double res;
    
    res=gggtt2virt(xp1,xp2,xp3,xkq,xkqb);
  
    //  cout<<RED<<res<<RESET<<endl;

    return res;
  }


   double ioperator_gg_(const double *p1,const double *p2, const double *p3,const double *kq,const double *kqb){
    FourMomentum xp1(p1), xp2(p2), xp3(p3), xkq(kq), xkqb(kqb);
    double res;

    res=Ioperator_gg(xp1,xp2,xp3,xkq,xkqb);
  
    // cout<<RED<<res<<RESET<<endl;

    return res;
  }





  // f77 interface  to Born in Uwer's library
  //
  // Corresponds to the following Fortran subroutine
  // gqbttqb2(p(0,1),p(0,2),p(0,5),p(0,3),p(0,4))
  // interface structure:
  //
  double gqbttqb2_(const double *p1,const double *p2, const double *p3,const double *kq,const double *kqb){
    FourMomentum xp1(p1), xp2(p2), xp3(p3), xkq(kq), xkqb(kqb);

    double res;
    
    res=gqbttqb2(xp1,xp2,xp3,xkq,xkqb);
  
    //  cout<<RED<<res<<RESET<<endl;

    return res;
  }

  // f77 interface to virtuals in Uwer's library
  //
  // Corresponds to the following Fortran subroutine
  // gqbttqb2virt(p(0,1),p(0,2),p(0,5),p(0,3),p(0,4))
  //
  double gqbttqb2virt_(const double *p1,const double *p2, const double *p3,const double *kq,const double *kqb){
    FourMomentum xp1(p1), xp2(p2), xp3(p3), xkq(kq), xkqb(kqb);

    double res;
    
    res=gqbttqb2virt(xp1,xp2,xp3,xkq,xkqb);
  
    //  cout<<RED<<res<<RESET<<endl;

    return res;
  }


   double ioperator_gqb_(const double *p1,const double *p2, const double *p3,const double *kq,const double *kqb){
    FourMomentum xp1(p1), xp2(p2), xp3(p3), xkq(kq), xkqb(kqb);
    double res;

    res=Ioperator_gqb(xp1,xp2,xp3,xkq,xkqb);
  
    // cout<<RED<<res<<RESET<<endl;

    return res;
  }



  // f77 interface to Born in Uwer's library
  //
  // Corresponds to the following Fortran subroutine
  // qgttq2(p(0,1),p(0,2),p(0,5),p(0,3),p(0,4))
  //
  double qgttq2_(const double *p1,const double *p2, const double *p3,const double *kq,const double *kqb){
    FourMomentum xp1(p1), xp2(p2), xp3(p3), xkq(kq), xkqb(kqb);

    double res;
    
    res=qgttq2(xp1,xp2,xp3,xkq,xkqb);
  
    ///cout<<RED<<res<<RESET<<endl;

    return res;
  }


  // f77 interface to virtuals in Uwer's library
  //
  // Corresponds to the following Fortran subroutine
  // qgttq2virt(p(0,1),p(0,2),p(0,5),p(0,3),p(0,4))
  //
  double qgttq2virt_(const double *p1,const double *p2, const double *p3,const double *kq,const double *kqb){
    FourMomentum xp1(p1), xp2(p2), xp3(p3), xkq(kq), xkqb(kqb);

    double res;
    
    res=qgttq2virt(xp1,xp2,xp3,xkq,xkqb);
  
    //  cout<<RED<<res<<RESET<<endl;

    return res;
  }


   double ioperator_qg_(const double *p1,const double *p2, const double *p3,const double *kq,const double *kqb){
    FourMomentum xp1(p1), xp2(p2), xp3(p3), xkq(kq), xkqb(kqb);
    double res;

    res=Ioperator_qg(xp1,xp2,xp3,xkq,xkqb);
  
    // cout<<RED<<res<<RESET<<endl;

    return res;
  }





  // f77 interface to Born in Uwer's library
  //
  // Corresponds to the following Fortran subroutine
  // qqttg2(p(0,1),p(0,2),p(0,5),p(0,3),p(0,4))
  //
  double qqttg2_(const double *p1,const double *p2, const double *p3,const double *kq,const double *kqb){
    FourMomentum xp1(p1), xp2(p2), xp3(p3), xkq(kq), xkqb(kqb);

    double res;
    
    res=qqttg2(xp1,xp2,xp3,xkq,xkqb);
  
    //  cout<<RED<<res<<RESET<<endl;

    return res;
  }


  // f77 interface to virtuals in Uwer's library
  //
  // Corresponds to the following Fortran subroutine
  // qqttg2virt(p(0,1),p(0,2),p(0,5),p(0,3),p(0,4))
  //
  double qqttg2virt_(const double *p1,const double *p2, const double *p3,const double *kq,const double *kqb){
    FourMomentum xp1(p1), xp2(p2), xp3(p3), xkq(kq), xkqb(kqb);

    double res;
    
    res=qqttg2virt(xp1,xp2,xp3,xkq,xkqb);
  
    //  cout<<RED<<res<<RESET<<endl;

    return res;
  }

   double ioperator_qq_(const double *p1,const double *p2, const double *p3,const double *kq,const double *kqb){
    FourMomentum xp1(p1), xp2(p2), xp3(p3), xkq(kq), xkqb(kqb);
    double res;

    res=Ioperator_qq(xp1,xp2,xp3,xkq,xkqb);
  
    // cout<<RED<<res<<RESET<<endl;

    return res;
  }

  
} //end of extern "C" 
