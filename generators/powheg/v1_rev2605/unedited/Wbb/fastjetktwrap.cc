//STARTHEADER
// $Id: fastjetsiscone.cc 939 2007-11-08 12:18:08Z salam $
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

#include <iostream>
#include "fastjet/ClusterSequence.hh"


namespace fj = fastjet;
using namespace std;

extern "C" {   

// f77 interface to fastjet written by SA&ER
//
// Corresponds to the following Fortran subroutine
// interface structure:
//
//   SUBROUTINE FASTJETKTWHICH(P,NPART,PTMIN,R,F77JETS,NJETS,F77JETVEC)
//   DOUBLE PRECISION P(4,*), PTMIN, R, F77JETS(4,*)
//   INTEGER          NPART, NJETS,F77JETVEC(*)
// 
// where on input
//
//   P           the input particle 4-momenta
//   NPART       the number of input momenta
//   PTMIN       the minimun pt of jets
//   R           the radius parameter
//
// and on output 
//
//   F77JETS     the output jet momenta (whose second dim should be >= NPART)
//               sorted in order of decreasing p_t.
//   NJETS       the number of output jets 
//   F77JETVEC   the jet in which the corresponding track was clustered

 void fastjetktwhich_(const double * p, const int & npart,                   
		      const double & ptmin, const double & R,                    
		      double * f77jets, int & njets, int * f77jetvec) {

    // transfer p[4*ipart+0..3] -> input_particles[i]
    vector<fj::PseudoJet> input_particles;   
    for (int i=0; i<npart; i++) {
      valarray<double> mom(4); // mom[0..3]
      for (int j=0;j<=3; j++) {
         mom[j] = *(p++);
      }
      input_particles.push_back(fj::PseudoJet(mom));
      // label input_particles entries
      input_particles[i].set_user_index(i+1);
      //cout<<"input particle index "<<input_particles[i].user_index()<<endl;
    }

    // useful to test the behaviour
    if(false) {
      cout<<"numb. of particles "<<input_particles.size()<<endl;
      for (int i=0; i<input_particles.size(); i++) {
	cout<<"track "<<i<<" p "<<input_particles[i].px()<<" "<<input_particles[i].py()<<" "
	    <<input_particles[i].pz()<<" "<<input_particles[i].E() /*<<" "<<input_particles[i].m()*/
	    <<endl;
      }
    }

    // prepare jet def and run fastjet
    fj::JetDefinition jet_def(fj::kt_algorithm, R);
    // perform clustering
    fj::ClusterSequence cs(input_particles,jet_def);
    // extract jets (pt-ordered)
    vector<fj::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptmin));
    //vector<fj::PseudoJet> jets = sorted_by_pt(cs.exclusive_jets(ptmin*ptmin));
    // extract number of jets
    njets = jets.size();
    //cout<<"njets in fastjet "<<njets<<endl;

    // find particles inside i-th jet
    vector<fj::PseudoJet> *constit;
    constit=new vector<fj::PseudoJet>[njets];
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

    delete [] constit;

   }
}

