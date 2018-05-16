// This wrapper is essentially the same as fastjetsisconewrap.cc

#include <iostream>
#include "fastjet/ClusterSequence.hh"
#include "fastjet/D0RunIIConePlugin.hh"

namespace fj = fastjet;
using namespace std;

extern "C" {   

// f77 interface to D0RunIICone (via fastjet)
//
// Corresponds to the following Fortran subroutine
// interface structure:
//
//   SUBROUTINE FASTJETD0RUNIICONE(P,NPART,R,ETMIN,F,F77JETS,NJETS)
//   DOUBLE PRECISION P(4,*), R, ETMIN, F, F77JETS(4,*)
//   INTEGER          NPART, NJETS
// 
// where on input
//
//   P        the input particle 4-momenta
//   NPART    the number of input momenta
//   R        the radius parameter
//   ETMIN    cones to be discarded at if at any iteration
//            they have pt < Et_min_ratio * min_jet_Et.
//   F        the overlap threshold
//
// and on output 
//
//   F77JETS  the output jet momenta (whose second dim should be >= NPART)
//            sorted in order of decreasing p_t.
//   NJETS    the number of output jets 
//
void fastjetd0runiicone_(const double * p, const int & npart,                   
                     const double & R, const double & ETMIN, const double & f,
			 double * f77jets, int & njets, double * f77_pT_rel_vec) {

    // transfer p[4*ipart+0..3] -> input_particles[i]
    vector<fj::PseudoJet> input_particles;   
    for (int i=0; i<npart; i++) {
      valarray<double> mom(4); // mom[0..3]
      for (int j=0;j<=3; j++) {
         mom[j] = *(p++);
      }
      fj::PseudoJet psjet(mom);
      input_particles.push_back(psjet);    
    }
    
    // prepare jet def and run fastjet
    fj::D0RunIIConePlugin * plugin = new fj::D0RunIIConePlugin(R,ETMIN,f);
    fj::JetDefinition jet_def(plugin);
    
    // perform clustering
    fj::ClusterSequence cs(input_particles,jet_def);
    // extract jets (pt-ordered)
    vector<fj::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
    njets = jets.size();

    // transfer jets -> f77jets[4*ijet+0..3]
    for (int i=0; i<njets; i++) {
      for (int j=0;j<=3; j++) {
        *f77jets = jets[i][j];
        f77jets++;
      } 
    }


    // Calculation of pT rel as defined in arXiv:0907.4076 (Alioli et al).

    for(unsigned int ixx=0; ixx<njets; ixx++) {
      // Calculate the mass of the longitudinal boost for the current jet:
      double zboost_vec_mass(jets[ixx].E ()*jets[ixx].E ()
			    -jets[ixx].pz()*jets[ixx].pz());
      if(zboost_vec_mass>0.) {
	zboost_vec_mass = sqrt(zboost_vec_mass);
	
	// Work out the parameters of the associated boost matrix:
	double zboost_gamma(jets[ixx].E()/zboost_vec_mass);
	double zboost_betagamma(jets[ixx].pz()/zboost_vec_mass);
	
	// Extract the constituents of the jet.
	vector<fj::PseudoJet> jet_constituents(cs.constituents(jets[ixx]));
	vector<fj::PseudoJet> jet_constituents_boosted;
	
	// Compute the pT of the jet.
	double jet_pT(sqrt(jets[ixx].kt2()));
	
	// Boost each of the particles in the jet to the same frame
	// and for each one compute pT rel:
	double pT_rel(0.);
	double current_pT_rel2(0.);
	fj::PseudoJet sum_constituents_boosted(0.,0.,0.,0.);
	for(unsigned int jxx=0; jxx<jet_constituents.size() ; jxx++ ) {
	  jet_constituents_boosted.push_back(
	    fj::PseudoJet(jet_constituents[jxx].px(),
			  jet_constituents[jxx].py(),
			  jet_constituents[jxx].pz()*zboost_gamma
			 -jet_constituents[jxx].E ()*zboost_betagamma, 
			  jet_constituents[jxx].E ()*zboost_gamma
			 -jet_constituents[jxx].pz()*zboost_betagamma
			 ));
	  current_pT_rel2 =
	       pow(jet_constituents_boosted[jxx].pz(),2)*pow(jet_pT,2)
	      +pow(jet_constituents_boosted[jxx].px()*jets[ixx].py()
	          -jet_constituents_boosted[jxx].py()*jets[ixx].px(),2);
	  if(current_pT_rel2>=0.) 
	    pT_rel += sqrt(current_pT_rel2);
	  else 
	    cout << "fastjetd0runiicone_\n"
		 << "Warning pT_rel^2 for particle " << jxx 
		 << " in jet " << ixx 
		 << " = " << current_pT_rel2 << endl
		 << "This particle is omitted in calculating the total "
		 << "pT_rel for this jet."
		 << endl;
	  
	  sum_constituents_boosted += jet_constituents_boosted[jxx];
	}
	pT_rel /= jet_pT;

        *f77_pT_rel_vec = pT_rel;
        f77_pT_rel_vec++;

	bool debugging(false);
	if(debugging) {
	  cout << "\n";
	  cout << "Before boost:\n";
	  cout << "jet = " 
	       << jets[ixx].px() << " " << jets[ixx].py() << " " 
	       << jets[ixx].pz() << " " << jets[ixx].E()  << " " << endl;
	  cout << "jet mass =  " << jets[ixx].m() << endl;
	  cout << "sum of boosted constituent momenta = " 
	       << sum_constituents_boosted.px() << " " 
	       << sum_constituents_boosted.py() << " " 
	       << sum_constituents_boosted.pz() << " " 
	       << sum_constituents_boosted.E()  << " " << endl;
	  cout << "sum of boosted constituent momenta mass = " 
	       << sum_constituents_boosted.m() << endl;
	  cout << "pT_rel = " << pT_rel << endl;
	}
      }
    }

    // clean up
    delete plugin;
    
   }
}

