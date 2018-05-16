#include <cstdio>
#include <iostream>
#include <cstdio>
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"


using namespace fastjet;
using namespace std;

static bool verbose=false;
//static bool verbose=true;

// a declaration of a function that pretty prints a list of jets
void print_jets (const ClusterSequence &,const vector<PseudoJet> &);
// a declaration the fastjet function, to be used inside the fortran wrapper
void fastjet_kt(const double *, const int &, const double &, 
		const double &, const double &, const int &, 
		const int &, const int &, 
		int &, double *, int *);


extern "C" {
// f77 interface to fastjet written by SA
//
// Corresponds to the following Fortran subroutine
// interface structure:
//
//   SUBROUTINE
//   FASTJET_KT(P,NPART,JETCUT,R,PALG,RECO,CLUSTER,SORT,F77JETS,NJETS,F77JETVEC)
//   DOUBLE PRECISION P(4,*), JETCUT, R, PALG, RECO,F77JETS(4,*) INTEGER
//   NPART, NJETS,F77JETVEC(*)
//   INTEGER CLUSTER,SORT
// 
// where on input
//
//   P           the input particle 4-momenta
//   NPART       the number of input momenta
//   JETCUT      the minimum pt of jets in Inclusive algo
//               the dcut value in Exclusive algo
//   R           the radius parameter
//   PALG        the power for the generalised kt alg 
//               (1.0=kt, 0.0=C/A,  -1.0 = anti-kt)
//   RECO        the recombination scheme
//               (0=E 1=pt 2=pt2 3=Et 4=Et2)
//   CLUSTER     the jet clustering algo  
//               (0=Inclusive JETCUT=ptmin
//                1=Exclusive JETCUT=dcut
//                2=Exclusive requiring just NJETS)
//   SORT        the sorting strategy (0=pt 1=E 2=y)
//  
// and on output 
//
//   NJETS       the number of output jets ( if CLUSTER=2 this is an input value) 
//   F77JETS     the output jet momenta (whose second dim should be >= NPART)
//   F77JETVEC   the jet in which the corresponding track was clustered
  void fastjet_kt_(const double * p, const int & npart,
		   const double & jetcut, const double & R, 
		   const double & Palg, const int & Reco, 
		   const int & Cluster, const int & Sort, 
		   int & njets, double * f77jets,int * f77jetvec) {   
    fastjet_kt(p,npart,jetcut,R,Palg,Reco,Cluster,Sort,njets,f77jets,f77jetvec) ;
  }
} // end of extern "C"

void fastjet_kt(const double * p, const int & npart, 
		   const double & jetcut, const double & R, 
		   const double & Palg, const int & Reco, 
		   const int & Cluster, const int & Sort, 
		   int & njets, double * f77jets,int * f77jetvec) {   
{

    // transfer p[4*ipart+0..3] -> input_particles[i]
    vector<PseudoJet> input_particles; 
    for (int i=0; i<npart; i++) {
      valarray<double> mom(4); // mom[0..3]
      for (int j=0;j<=3; j++) {
         mom[j] = *(p++);
      }
      input_particles.push_back(PseudoJet(mom));
      // label input_particles entries
      input_particles[i].set_user_index(i+1);
      //cout<<"input particle index "<<input_particles[i].user_index()<<endl;
    }

    // useful to test the behaviour
    if(verbose) {
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

    //Select recombination strategy
    RecombinationScheme recomb_scheme;
    switch(Reco) {
    case 0:
      recomb_scheme = E_scheme;
      break;
    case 1:
      recomb_scheme = pt_scheme;
      break;
    case 2:
      recomb_scheme = pt2_scheme;
      break;
    case 3:
      recomb_scheme = Et_scheme;
      break;
    case 4:
      recomb_scheme = Et2_scheme;
      break;
    default:
      cout<<"Recombination scheme not allowed yet"<<endl;
      exit(1);
    }

     
    // prepare jet def and run fastjet
    JetDefinition jet_def;
    if (Palg == 1.0) {
      jet_def = JetDefinition(kt_algorithm, R,recomb_scheme, strategy);
    }  else if (Palg == 0.0) {
      jet_def = JetDefinition(cambridge_algorithm, R,recomb_scheme, strategy);
    }  else if (Palg == -1.0) {
      jet_def = JetDefinition(antikt_algorithm, R,recomb_scheme, strategy);
    } else {
      jet_def = JetDefinition(genkt_algorithm, R, Palg,recomb_scheme, strategy);
    }

    
    // perform clustering
    ClusterSequence cs(input_particles, jet_def);

    // extract jets 
    vector<PseudoJet> jets;
    if (Cluster == 0) {
      jets= cs.inclusive_jets(jetcut);
      njets = jets.size();
    } else if (Cluster == 1) {
      jets= cs.exclusive_jets(jetcut);
      njets = jets.size();
    } else if (Cluster == 2) {
      jets= cs.exclusive_jets(njets);
    }
    
    // sort
    switch (Sort) {
    case 0:
      jets=sorted_by_pt(jets);
      break;
    case 1:
      jets=sorted_by_E(jets);
      break;
    case 2: 
      jets=sorted_by_rapidity(jets);
      break;
    }
    
    static bool first = true; 
    if(first){
      cout << "# -------------------------------------------------------------" 
	   << endl;
      cout << "# Running " << jet_def.description() << endl;
      cout << "# Strategy adopted by FastJet is "<<
	cs.strategy_string()<<endl;
      if(Cluster==0) {
	cout << "# Accessing Inclusive Jets  ptmin="<<jetcut<<endl; }
      else if(Cluster==1) {
	cout << "# Accessing Exclusive Jets  dcut="<<jetcut<<endl; }         
      else { 
	cout << "# Accessing Exclusive Jets  njets="<<njets<<endl; }     
      if(Sort==3) {
	cout << "# Jets sorted by rapidity"<<endl; }
      else if (Sort==2){
	cout << "# Jets sorted by E"<<endl; }         
      else { 
	cout << "# Jets sorted by pt"<<endl; }   
      first=false;
  cout << "# -------------------------------------------------------------" 
       << endl;
  cout << "#" << endl;
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

 if(verbose) {
      // print them out
      cout << "Printing jets \n";
      cout << "---------------------------------------\n";
      print_jets(cs, jets);
      cout << endl;
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
}


/// a function that pretty prints a list of jets
void print_jets (const ClusterSequence & clust_seq, 
		 const vector<PseudoJet> & jets) {

  // label the columns
  printf("%5s %15s %15s %15s %15s\n","jet #", "rapidity", 
	 "phi", "pt", "n constituents");
  
  // print out the details for each jet
  for (unsigned int i = 0; i < jets.size(); i++) {
    int n_constituents = clust_seq.constituents(jets[i]).size();
    printf("%5u %15.8f %15.8f %15.8f %8u\n",
	   i, jets[i].rap(), jets[i].phi(),
	   jets[i].perp(), n_constituents);
  }

}



