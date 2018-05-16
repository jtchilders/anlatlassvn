#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/MassDropTagger.hh"

namespace fj = fastjet;
using namespace fastjet;
using namespace std;



extern "C" {  
  
  void fastjetfatjet_(const double * p, const int & npart,                   
		      const double & R, const double & palg, const double & ptmin,
		      double * f77jets, int & njets, int * f77jetvec, int & itagjet, double * f77mwjet, 
		      bool & passcuts) {
    
    // transfer p[4*ipart+0..3] -> input_particles[i]
    vector<fj::PseudoJet> input_particles;   
    for (int i=0; i<npart; i++) {
      valarray<double> mom(4); // mom[0..3]
      for (int j=0;j<=3; j++) {
	mom[j] = *(p++);
      }
      fj::PseudoJet psjet(mom);
      input_particles.push_back(psjet);    
      // label input_particles entries
      input_particles[i].set_user_index(i+1);
    }
      
    // prepare jet def and run fastjet
    fj::JetDefinition jet_def;
    if (palg == 1.0) {
      jet_def = fj::JetDefinition(fj::kt_algorithm, R);
    }  else if (palg == 0.0) {
      jet_def = fj::JetDefinition(fj::cambridge_algorithm, R);
    }  else if (palg == -1.0) {
      jet_def = fj::JetDefinition(fj::antikt_algorithm, R);
    } else {
      jet_def = fj::JetDefinition(fj::genkt_algorithm, R, palg);
    }
    
    
    // perform clustering
    fj::ClusterSequence cs(input_particles, jet_def);
    // extract jets (pt-ordered)
    vector<fj::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptmin));
    njets = jets.size();

    // initialize to 0 
    itagjet = -1; 
    passcuts = true ; 

    // find particles inside i-th jet
    vector<fj::PseudoJet> *constit;
    constit=new vector<fj::PseudoJet>[njets];
    for(int i=0; i<njets; i++) {
      constit[i] = cs.constituents(jets[i]); 
      for(int j=0; j<constit[i].size(); j++) {
	*(f77jetvec + constit[i][j].user_index()-1) = i+1;
      }
    }

    // GZ -- NEW STUFF AFTER HERE 
    double mw = 80.419; 
    double mmin = 1000000000; 
    PseudoJet fatjet; 
    double masswindow = 10; 
    double ptcut2 = 300*300; 
    double ptmin2 = 0; 
    // declare the tagger 
    double mcut = 0.67; 
    double ycut = 0.09; 
    MassDropTagger mdtagger(mcut, ycut);

    vector<fj::PseudoJet> tagged_jets = mdtagger(jets);
    //    PseudoJet filtered_tagged_jet;
    for(int i=0; i<njets; i++) {
      if (tagged_jets[i] != 0){
	
	// output filtered jet 
	vector<PseudoJet> pieces = tagged_jets[i].pieces(); 
	double Rfilt = min(0.3, 0.5 * pieces[0].delta_R(pieces[1])); 
	PseudoJet filtered_tagged_jet = Filter(Rfilt, SelectorNHardest(3))(tagged_jets[i]);
	
//	if (abs(filtered_tagged_jet.m2()-mw2) <  mmin && filtered_tagged_jet.perp2() > ptcut2 ){
//	  mmin  = abs(filtered_tagged_jet.m2()-mw2); 
//	  fatjet = filtered_tagged_jet; 
//	  itagjet = i; 
//	}

	if (abs(sqrt(filtered_tagged_jet.m2())-mw) <  masswindow  
	    && filtered_tagged_jet.perp2() > ptmin2 
	    && filtered_tagged_jet.perp2() > ptcut2){
	  ptmin2 = filtered_tagged_jet.perp2(); 
	  fatjet = filtered_tagged_jet; 
	  itagjet = i; 
	}

      }
    }
    
    if (itagjet == -1){
      passcuts = false ;
    }

//    // now require jet to be in the 10 GeV mass window 
//    if( abs(fatjet.m2()-mw2) > masswindow2) {
//      passcuts = false ; 
//    }

    if (passcuts) {
      // transfer jets -> f77jets[4*ijet+0..3]
      for (int i=0; i<njets; i++) {
	for (int j=0;j<=3; j++) {
	  *f77jets = jets[i][j];
	  f77jets++;
	} 
      }
      
      for (int j=0;j<=3; j++) {
	*f77mwjet = fatjet[j];
	//cout << "filt" << filtered_tagged_jet[j] << endl; 
	*f77mwjet++;
      }
      // to fortran convention (1...n) rather than (0...n-1) 
      itagjet=itagjet+1;

    } // if passcuts
    //    cout << "passcuts" << passcuts << endl; 

      // END OF NEW STUFF 
      // clean up
    delete [] constit;
    
  }
}

