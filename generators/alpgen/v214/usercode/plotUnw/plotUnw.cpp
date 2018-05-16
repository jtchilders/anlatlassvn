// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk. 
// It studies the charged multiplicity distribution at the LHC.

#include "GeneratorInput.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"
#include "tclap/ArgException.h"

double get_pT(double px,double py){
   return sqrt(px*px + py*py);
}

std::vector<unsigned int>* getPtOrderedIndices(LHAupAlpgen* lha,std::vector<unsigned int>* unordered_indices){

   std::vector<unsigned int>* ordered_indices = new std::vector<unsigned int>;


   for(std::vector<unsigned int>::iterator it = unordered_indices->begin();it != unordered_indices->end();++it){
      double pT = get_pT(lha->px(*it),lha->py(*it));

      // std::cout << " > " << (it-unordered_indices->begin()) << " " << *it << " " << pT << "\n";

      // first element can be just pushed back
      if(ordered_indices->size() == 0){
         ordered_indices->push_back(*it);
      }
      else{
         // loop over ordered indices to find where this element should be inserted
         unsigned int initial_size = ordered_indices->size();
         for(std::vector<unsigned int>::iterator ordered_iter = ordered_indices->begin();
              ordered_iter != ordered_indices->end();
              ++ordered_iter){
            double ordered_pT = get_pT(lha->px(*ordered_iter),lha->py(*ordered_iter));

            if(pT > ordered_pT){
               ordered_indices->insert(ordered_iter,*it);
               break;
            }
         }
         // if vector hasn't grown, push back value
         if(initial_size == ordered_indices->size()){
            ordered_indices->push_back(*it);
         }

      }
   }

   // for(unsigned int i=0;i<ordered_indices->size();++i)
   //    std::cout << " < " << i << " " << ordered_indices->at(i) << "\n";

   return ordered_indices;
}

class ParticlePlots{
public:
   ParticlePlots(std::string particle_name){
      m_particle_name = particle_name;
      std::stringstream name,title;
      double X_MIN = 0.;
      double X_MAX = 3000.;
      double X_BINS = 100;
      
      name << "th1_" << m_particle_name << "_px";
      title << ";" << m_particle_name << " p_{x} [GeV]; events";
      m_h_px = new TH1D(name.str().c_str(),title.str().c_str(),X_BINS,X_MIN,X_MAX);

      name.str("");title.str("");
      name << "th1_" << m_particle_name << "_py";
      title << ";" << m_particle_name << " p_{y} [GeV]; events";
      m_h_py = new TH1D(name.str().c_str(),title.str().c_str(),X_BINS,X_MIN,X_MAX);

      name.str("");title.str("");
      name << "th1_" << m_particle_name << "_pz";
      title << ";" << m_particle_name << " p_{z} [GeV]; events";
      m_h_pz = new TH1D(name.str().c_str(),title.str().c_str(),X_BINS,X_MIN,X_MAX);

      name.str("");title.str("");
      name << "th1_" << m_particle_name << "_energy";
      title << ";" << m_particle_name << " Energy [GeV]; events";
      m_h_energy = new TH1D(name.str().c_str(),title.str().c_str(),X_BINS,X_MIN,X_MAX);

      name.str("");title.str("");
      name << "th1_" << m_particle_name << "_mass";
      title << ";" << m_particle_name << " Mass [GeV]; events";
      m_h_mass = new TH1D(name.str().c_str(),title.str().c_str(),X_BINS,X_MIN,X_MAX);

   }
   ~ParticlePlots(){
      delete m_h_px; m_h_px = 0;
      delete m_h_py; m_h_py = 0;
      delete m_h_pz; m_h_pz = 0;
      delete m_h_energy; m_h_energy = 0;
      delete m_h_mass; m_h_mass = 0;
   }

   void fillParticleHistograms(LHAupAlpgen const * const lha,unsigned int index){
      if(index >= lha->sizePart()){
         std::cerr << " Error: Index passed, " << index << ", is greater than the number of particles in the array, " << lha->sizePart() << "\n";
         return;
      }
      m_h_px->Fill(lha->px(index));
      m_h_py->Fill(lha->py(index));
      m_h_pz->Fill(lha->pz(index));
      m_h_energy->Fill(lha->e(index));
      m_h_mass->Fill(lha->m(index));
   }

   void Write(){
      m_h_px->Write(m_h_px->GetName(),TObject::kOverwrite);
      m_h_py->Write(m_h_py->GetName(),TObject::kOverwrite);
      m_h_pz->Write(m_h_pz->GetName(),TObject::kOverwrite);
      m_h_energy->Write(m_h_energy->GetName(),TObject::kOverwrite);
      m_h_mass->Write(m_h_mass->GetName(),TObject::kOverwrite);
   }
private:
   TH1D* m_h_px;
   TH1D* m_h_py;
   TH1D* m_h_pz;
   TH1D* m_h_energy;
   TH1D* m_h_mass;
   std::string m_particle_name;
};

std::vector<unsigned int>* getParticleIdIndex(LHAupAlpgen* const lha,int particle_id,int status_req = 1){
   std::vector<unsigned int>* particles = new std::vector<unsigned int>;
   for(unsigned int i = 0;i<lha->sizePart();++i){
      if( lha->id(i) == particle_id && lha->status(i) == status_req)
         particles->push_back(i);
   }
   return particles;
}

std::vector<unsigned int>* getAbsParticleIdIndex(LHAupAlpgen* const lha,int particle_id,int status_req = 1){
   std::vector<unsigned int>* particles = new std::vector<unsigned int>;
   for(unsigned int i = 0;i<lha->sizePart();++i){
      if( fabs(lha->id(i)) == fabs(particle_id) && lha->status(i) == status_req)
         particles->push_back(i);
   }
   return particles;
}

std::vector<unsigned int>* getParticleIdIndex(LHAupAlpgen* const lha,std::vector<int>& particle_ids,int status_req = 1){
   std::vector<unsigned int>* particles = new std::vector<unsigned int>;
   for(unsigned int i = 0;i<lha->sizePart();++i){
      const int pid = lha->id(i);
      std::vector<int>::const_iterator it = std::find(particle_ids.begin(),particle_ids.end(),pid);
      if( it != particle_ids.end() && lha->status(i) == status_req)
         particles->push_back(i);
   }
   return particles;
}

std::vector<unsigned int>* getAbsParticleIdIndex(LHAupAlpgen* const lha,std::vector<int>& particle_ids,int status_req = 1){
   std::vector<unsigned int>* particles = new std::vector<unsigned int>;
   for(unsigned int i = 0;i<lha->sizePart();++i){
      const int pid = lha->id(i);
      std::vector<int>::const_iterator it = std::find(particle_ids.begin(),particle_ids.end(),fabs(pid));
      if( it != particle_ids.end() && lha->status(i) == status_req)
         particles->push_back(i);
      it = std::find(particle_ids.begin(),particle_ids.end(),-1*fabs(pid));
      if( it != particle_ids.end() && lha->status(i) == status_req)
         particles->push_back(i);
   }
   return particles;
}


int main(int argc,char** argv) {

   std::string input_filename,output_filename;
   bool quiet = false;
   try{
      TCLAP::CmdLine cmd(" Join multiple HepMC files into one ",' ',"0.1");
      TCLAP::ValueArg<std::string> alpgenBaseArg("i","alpgen-base","Alpgen basename such as, 'alpout', which will be used to find the .unw and _unw.par files to plot.",true,"","string",cmd);
      TCLAP::ValueArg<std::string> outputFileArg("o","outfilename","Output filename",false,"output.root","string",cmd);
      TCLAP::SwitchArg quietArg("q","quiet","Do not count events",cmd);
      cmd.parse(argc,argv);
      input_filename = alpgenBaseArg.getValue();
      output_filename = outputFileArg.getValue();
      quiet = quietArg.getValue();
   }
   catch (TCLAP::ArgException &e){
      std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
   }
   
   std::cout << " alpgen base filename: " << input_filename << "\n";
   std::cout << " output filename: " << output_filename << "\n";
   
   // jet PIDs
   std::vector<int> jet_pids;
   // jets only, so quarks(1 <= pdg id <= 6) and gluons(pdg id == 21)
   jet_pids.push_back(1); // d
   jet_pids.push_back(2); // u 
   jet_pids.push_back(3); // s
   jet_pids.push_back(4); // c
   jet_pids.push_back(5); // b
   jet_pids.push_back(6); // t
   jet_pids.push_back(21); // gluon

   // open alpgen input files
   LHAupAlpgen* lha = new LHAupAlpgen(input_filename.c_str());
   //std::string output_filename = "output.root";
   //LHAupAlpgen* lha = new LHAupAlpgen("alpout");
   if(!lha->fileFound()){ std::cerr << " File Not Found \n";return -1;}

   // initialize alpgen event reader
   lha->setInit();

   // initialize plots
   std::vector<ParticlePlots*> pp_jets;
   ParticlePlots* pp_z = 0;
   ParticlePlots* pp_w = 0;
   ParticlePlots* pp_el_plus = 0;
   ParticlePlots* pp_el_minus = 0;
   ParticlePlots* pp_mu_plus = 0;
   ParticlePlots* pp_mu_minus = 0;
   ParticlePlots* pp_tau_plus = 0;
   ParticlePlots* pp_tau_minus = 0;
   ParticlePlots* pp_nu_el = 0;
   ParticlePlots* pp_nu_mu = 0;
   ParticlePlots* pp_nu_tau = 0;

   
   // loop over events (setEvent function returns true until end of event list reached)
   int i = 0;
   int max_i = -1; // can limit the number of events looped over.
   while(lha->setEvent(0,-1) && (i++ < max_i || max_i == -1)) {
      if( i % 10000 == 0 and !quiet) std::cout << " Processed  " << i << " Events.\n";

      // for(unsigned int p=0;p<lha->sizePart();++p){
      //    std::cerr << " > " << p << " pid: " << lha->id(p) << " status: " << lha->status(p) << " E: " << lha->e(p) << " M: " << lha->m(p) << "\n";
      // }

      // plot any jets, using pT ordering of the jets
      // first get a list of the vector indices within the lhaparticle vector for this event
      std::vector<unsigned int>* unordered_jet_particle_indices = getAbsParticleIdIndex(lha,jet_pids);
      std::vector<unsigned int>* ordered_jet_particle_indices = getPtOrderedIndices(lha,unordered_jet_particle_indices);
      for(unsigned int i=0;i<ordered_jet_particle_indices->size();++i){
         unsigned int index = ordered_jet_particle_indices->at(i);
         
         // do I need to add a plot to the pp_jets
         ParticlePlots* pp = 0;
         if(pp_jets.size() <= i){
            std::stringstream name;
            name << "jet_" << i;
            pp = new ParticlePlots(name.str());
            pp_jets.push_back(pp);
         }
         else pp = pp_jets[i];
         //std::cerr << " i = " << i << " particle index = " << index << " pt = " << v->Pt() << "\n"; 
         
         // Fill plots
         pp->fillParticleHistograms(lha,index);
      }

      // plot any Z-bosons
      std::vector<unsigned int>* z_boson_index = getAbsParticleIdIndex(lha,23,2);
      for(std::vector<unsigned int>::iterator it = z_boson_index->begin();it != z_boson_index->end();++it){
         if(pp_z == 0){
            std::stringstream name;
            name << "z_boson";
            pp_z = new ParticlePlots(name.str());
         }
         pp_z->fillParticleHistograms(lha,*it);
      }

      // plot any W-bosons
      std::vector<unsigned int>* w_boson_index = getAbsParticleIdIndex(lha,24,2);
      for(std::vector<unsigned int>::iterator it = w_boson_index->begin();it != w_boson_index->end();++it){
         if(pp_w == 0){
            std::stringstream name;
            name << "w_boson";
            pp_w = new ParticlePlots(name.str());
         }
         pp_w->fillParticleHistograms(lha,*it);
      }

      // plot any e-
      std::vector<unsigned int>* index = getParticleIdIndex(lha,11);
      for(std::vector<unsigned int>::iterator it = index->begin();it != index->end();++it){
         if(pp_el_minus == 0){
            std::stringstream name;
            name << "e_minus";
            pp_el_minus = new ParticlePlots(name.str());
         }
         pp_el_minus->fillParticleHistograms(lha,*it);
      }

      // plot any e+
      index = getParticleIdIndex(lha,-11);
      for(std::vector<unsigned int>::iterator it = index->begin();it != index->end();++it){
         if(pp_el_plus == 0){
            std::stringstream name;
            name << "e_plus";
            pp_el_plus = new ParticlePlots(name.str());
         }
         pp_el_plus->fillParticleHistograms(lha,*it);
      }

      // plot any electon neutrinos
      index = getAbsParticleIdIndex(lha,12);
      for(std::vector<unsigned int>::iterator it = index->begin();it != index->end();++it){
         if(pp_nu_el == 0){
            std::stringstream name;
            name << "electron_neutrino";
            pp_nu_el = new ParticlePlots(name.str());
         }
         pp_nu_el->fillParticleHistograms(lha,*it);
      }

      // plot any mu-
      index = getParticleIdIndex(lha,13);
      for(std::vector<unsigned int>::iterator it = index->begin();it != index->end();++it){
         if(pp_mu_minus == 0){
            std::stringstream name;
            name << "mu_minus";
            pp_mu_minus = new ParticlePlots(name.str());
         }
         pp_mu_minus->fillParticleHistograms(lha,*it);
      }

      // plot any mu+
      index = getParticleIdIndex(lha,-13);
      for(std::vector<unsigned int>::iterator it = index->begin();it != index->end();++it){
         if(pp_mu_plus == 0){
            std::stringstream name;
            name << "mu_plus";
            pp_mu_plus = new ParticlePlots(name.str());
         }
         pp_mu_plus->fillParticleHistograms(lha,*it);
      }

      // plot any muon neutrinos
      index = getAbsParticleIdIndex(lha,12);
      for(std::vector<unsigned int>::iterator it = index->begin();it != index->end();++it){
         if(pp_nu_mu == 0){
            std::stringstream name;
            name << "muon_neutrino";
            pp_nu_mu = new ParticlePlots(name.str());
         }
         pp_nu_mu->fillParticleHistograms(lha,*it);
      }

      // plot any tau-
      index = getParticleIdIndex(lha,15);
      for(std::vector<unsigned int>::iterator it = index->begin();it != index->end();++it){
         if(pp_tau_minus == 0){
            std::stringstream name;
            name << "tau_minus";
            pp_tau_minus = new ParticlePlots(name.str());
         }
         pp_tau_minus->fillParticleHistograms(lha,*it);
      }

      // plot any tau+
      index = getParticleIdIndex(lha,-15);
      for(std::vector<unsigned int>::iterator it = index->begin();it != index->end();++it){
         if(pp_tau_plus == 0){
            std::stringstream name;
            name << "tau_plus";
            pp_tau_plus = new ParticlePlots(name.str());
         }
         pp_tau_plus->fillParticleHistograms(lha,*it);
      }

      // plot any tau neutrinos
      index = getAbsParticleIdIndex(lha,12);
      for(std::vector<unsigned int>::iterator it = index->begin();it != index->end();++it){
         if(pp_nu_tau == 0){
            std::stringstream name;
            name << "tau_neutrino";
            pp_nu_tau = new ParticlePlots(name.str());
         }
         pp_nu_tau->fillParticleHistograms(lha,*it);
      }

      
   }

   TFile* file = TFile::Open(output_filename.c_str(),"RECREATE");
   
   // write jet plots to file
   for(unsigned int i=0;i<pp_jets.size();++i)
      if(pp_jets[i]) pp_jets[i]->Write();
   
   // write others
   if(pp_z) pp_z->Write();
   if(pp_w) pp_w->Write();
   if(pp_el_minus) pp_el_minus->Write();
   if(pp_el_plus) pp_el_plus->Write();
   if(pp_nu_el) pp_nu_el->Write();
   if(pp_mu_minus) pp_mu_minus->Write();
   if(pp_mu_plus) pp_mu_plus->Write();
   if(pp_nu_mu) pp_nu_mu->Write();
   if(pp_tau_minus) pp_tau_minus->Write();
   if(pp_tau_plus) pp_tau_plus->Write();
   if(pp_nu_tau) pp_nu_tau->Write();

   file->Close();

   return 0;
}




