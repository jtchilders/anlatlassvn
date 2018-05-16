#include <iostream>
#include "fastjet/ClusterSequence.hh"
#include "fastjet/CDFMidPointPlugin.hh"

namespace fj = fastjet;
using namespace std;

extern "C" {   

// f77 interface to CDF MidPoint algorithm (via fastjet)
//
// Corresponds to the following Fortran subroutine
// interface structure:
//
//   SUBROUTINE FASTJETCDFMIDPOINT(P,NPART,R,F,SF,CAF,F77JETS,NJETS)
//   DOUBLE PRECISION P(4,*), R,F,SF,CAF, F77JETS(4,*)
//   INTEGER          NPART, NJETS
// 
// where on input
//
//   P        the input particle 4-momenta
//   NPART    the number of input momenta
//   R        the radius parameter
//   F        the overlap threshold, usually 0.5 in CDF 
//   SF       the seed threshold =1 in CDF 
//   CAF      the cone area fraction =1 in CDF
//
// and on output 
//
//   F77JETS  the output jet momenta (whose second dim should be >= NPART)
//            sorted in order of decreasing p_t.
//   NJETS    the number of output jets 
//
void fastjetcdfmidpoint_(const double * p, const int & npart,                   
                     const double & R, const double & f, const double & sf,                  
                     const double & caf, double * f77jets, int & njets) {

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
    fj::CDFMidPointPlugin * plugin = new fj::CDFMidPointPlugin(R,f,sf,caf);
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

    // clean up
    delete plugin;
    
   }
}

