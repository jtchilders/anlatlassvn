//STARTHEADER
// $Id: fastjetfortran.cc 1570 2009-05-25 10:45:18Z salam $
//
// Copyright (c) 2005-2007, Matteo Cacciari, Gavin Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet; if not, write to the Free Software
//  Foundation, Inc.:
//      59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//----------------------------------------------------------------------
//ENDHEADER
#include <cstdio>
#include <iostream>
#include <vector> 
#include <memory>
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"

using namespace fastjet;
using namespace std;

// a declaration of a function that pretty prints a list of jets
void print_jets (const ClusterSequence &,const vector<PseudoJet> &);

// a declaration the fastjet function, to be used inside the fortran wrapper
void fastjetppgenkt(const double * , const int & ,                   
		    const double &,  const double & , const double &,
		     double * , int &,  int * ) ;

extern "C" {   
// f77 interface to the pp generalised-kt (sequential recombination)
// algorithms, as defined in arXiv.org:0802.1189, which includes
// kt, Cambridge/Aachen and anti-kt as special cases.
//
// Corresponds to the following Fortran subroutine
// interface structure:
//
//   SUBROUTINE FASTJETPPSEQREC(P,NPART,R,PALG,F77JETS,NJETS)
//   DOUBLE PRECISION P(4,*), R, PALG, F, F77JETS(4,*)
//   INTEGER          NPART, NJETS
// 
// where on input
//
//   P        the input particle 4-momenta
//   NPART    the number of input momenta
//   R        the radius parameter
//   PALG     the power for the generalised kt alg 
//            (1.0=kt, 0.0=C/A,  -1.0 = anti-kt)
//
// and on output 
//
//   F77JETS  the output jet momenta (whose second dim should be >= NPART)
//            sorted in order of decreasing p_t.
//   NJETS    the number of output jets 
//
// For the values of PALG that correspond to "standard" cases (1.0=kt,
// 0.0=C/A, -1.0 = anti-kt) this routine actually calls the direct
// implementation of those algorithms, whereas for other values of
// PALG it calls the generalised kt implementation.
//
void fastjetppgenkt_(const double * p, const int & npart,                   
                     const double & R, const double & palg, const double & ptmin,
                     double * f77jets, int & njets, int * f77jetvec) {
  fastjetppgenkt(p,npart,R,palg,ptmin,f77jets,njets,f77jetvec);
}
}

void fastjetppgenkt(const double * p, const int & npart,                   
                     const double & R, const double & palg, const double & ptmin,
                     double * f77jets, int & njets, int * f77jetvec) {
    // transfer p[4*ipart+0..3] -> input_particles[i]
    vector<PseudoJet> input_particles;   
    for (int i=0; i<npart; i++) {
      valarray<double> mom(4); // mom[0..3]
      for (int j=0;j<=3; j++) {
         mom[j] = *(p++);
      }
      PseudoJet psjet(mom);
      input_particles.push_back(psjet);    
      // label input_particles entries
      input_particles[i].set_user_index(i+1);
    }
    
    // useful to test the behaviour
    if(false) {
      cout<<"numb. of particles "<<input_particles.size()<<endl;
      for (int i=0; i<input_particles.size(); i++) {
	cout<<"track "<<i<<" p "<<input_particles[i].px()<<" "<<input_particles[i].py()<<" "
	    <<input_particles[i].pz()<<" "<<input_particles[i].E() <<" "<<input_particles[i].m()
	    <<endl;
      }
    }

     // create an object that represents your choice of jet algorithm and 
    // the associated parameters
    double Rparam = R;
    Strategy strategy = Best;
    

    // prepare jet def and run fastjet
    JetDefinition jet_def;
    if (palg == 1.0) {
      jet_def = JetDefinition(kt_algorithm, R);
    }  else if (palg == 0.0) {
      jet_def = JetDefinition(cambridge_algorithm, R);
    }  else if (palg == -1.0) {
      jet_def = JetDefinition(antikt_algorithm, R);
    } else {
      jet_def = JetDefinition(genkt_algorithm, R, palg);
    }

    // perform clustering
    ClusterSequence cs(input_particles, jet_def);
     static bool first = true; 
    if(first){
      cout << "# -------------------------------------------------------------" 
       << endl;
      cout << "# Running " << jet_def.description() << endl;
      cout << "# Strategy adopted by FastJet is "<<
	cs.strategy_string()<<endl;
      first=false;
  cout << "# -------------------------------------------------------------" 
       << endl;
  cout << "#" << endl;
    }


    // extract jets (pt-ordered)
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptmin));
    njets = jets.size();

    if(false) {
      // print them out
      cout << "Printing inclusive jets with pt > "<< ptmin<<" GeV\n";
      cout << "---------------------------------------\n";
      print_jets(cs, jets);
      cout << endl;
    }


    // find particles inside i-th jet
    vector<PseudoJet> *constit;
    constit=new vector<PseudoJet>[njets];
    for(int i=0; i<njets; i++) {
      constit[i] = cs.constituents(jets[i]); 
      //cout<<"jet "<<i<<endl;
      //cout<<"mult "<<constit[i].size()<<endl;
      for(int j=0; j<constit[i].size(); j++) {
	*(f77jetvec + constit[i][j].user_index()-1) = i+1;
      }
    }



    // transfer jets -> f77jets[4*ijet+0..3]
    for (int i=0; i<njets; i++) {
      for (int j=0;j<=3; j++) {
        *f77jets = jets[i][j];
        f77jets++;
      } 
    }

    // clean up
    delete [] constit;
    
}

/// a function that pretty prints a list of jets
void print_jets (const ClusterSequence & clust_seq, 
		 const vector<PseudoJet> & jets) {

  // sort jets into increasing pt
  vector<PseudoJet> sorted_jets = sorted_by_pt(jets);  

  // label the columns
  printf("%5s %15s %15s %15s %15s\n","jet #", "rapidity", 
	 "phi", "pt", "n constituents");
  
  // print out the details for each jet
  for (unsigned int i = 0; i < sorted_jets.size(); i++) {
    int n_constituents = clust_seq.constituents(sorted_jets[i]).size();
    printf("%5u %15.8f %15.8f %15.8f %8u\n",
	   i, sorted_jets[i].rap(), sorted_jets[i].phi(),
	   sorted_jets[i].perp(), n_constituents);
  }

}

