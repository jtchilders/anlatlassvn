///////////////////////////////////////////////
//// Shower Alpgen output with Pythai 6
//// and save output as HepMC format  
//// written by Taylor Childers
//// based on pyuser example from Alpgen(2.14) 
//// and HepMC/example/fio/example_MyPythia.cc
///////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>

#include "HepMC/PythiaWrapper.h"
#include "HepMC/IO_HEPEVT.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/IO_AsciiParticles.h"
#include "HepMC/GenEvent.h"

int get_num_events(std::string filename);

int main(int argc,char** argv){
   //........................................HEPEVT
   // Pythia 6.1 uses HEPEVT with 4000 entries and 8-byte floating point
   //  numbers. We need to explicitly pass this information to the
   //  HEPEVT_Wrapper.
   //
   HepMC::HEPEVT_Wrapper::set_max_number_entries(4000);
   HepMC::HEPEVT_Wrapper::set_sizeof_real(8);
  
   ///////////////////////////////
   // initialize pythia
   //////////////////////

   // set Perugia Tune
   pypars.mstp[5-1]        = 356;
   // other ATLAS parameters
   pydat2.pmas[6-1][1-1]   = 172.5;    // top mass
   pydat2.pmas[23-1][1-1]  = 91.1876;  // Z mass
   pydat2.pmas[24-1][1-1]  = 80.399;   // W mass
   pydat1.paru[102-1]      = 0.23133;  // Weak Mixing Angle
   pydat1.parj[90-1]       = 20000.;
   pydat3.mdcy[15-1][1-1]  = 0;
   pypars.mstp[143-1]      = 1;        // UPVETO can be called

   call_pyinit("USER"," "," ",0.);
   
   // Instantiate an IO strategy for reading from HEPEVT.
   HepMC::IO_HEPEVT hepevtio;
   
   // Instantiate an IO strategy to write the data to file
   HepMC::IO_GenEvent ascii_io("pythia_output.hepmc",std::ios::out);
  
   // event loop
   std::cout << "\n";
   for(unsigned evt_num = 0;evt_num<eventMax;++evt_num){
      if(evt_num%50==1){
         std::stringstream s;
         s << "\n Event number: ";
         s.width(10);
         s << evt_num << " of ";
         s.width(10);
         s << eventMax << '\n';

         std::cout << s.str();
      }
      // shower event
      call_pyevnt();
      // convert from pythia PYJET format to HEPEVT
      call_pyhepc(1);
      // convert HEPEVT to HemMc format
      HepMC::GenEvent* evt = hepevtio.read_next_event();
      // need to set more parameters
      evt->set_event_number(evt_num);
      evt->set_signal_process_id(pysubs.msub[20-1]);
      // set number of multi parton interactions
      evt->set_mpi( pypars.msti[31-1] );
      // set cross section
      evt->set_cross_section( HepMC::getPythiaCrossSection() );
      // write the event out to the ascii files
      ascii_io << evt;
      // we also need to delete the created event from memory
      delete evt;
   }
   std::cout << "\n";
   call_pystat(1);
   
   return 0;
}

int get_num_events(std::string filename){
   
   std::fstream fin(filename.c_str(),std::fstream::in);

   std::string line;
   int numberEvents = -1;
   while(std::getline(fin,line)){
      if(line.find("! unwtd events, lum (pb-1)") >= 0){
         std::stringstream s(line);
         s >> numberEvents;
      }
      
   }
   return numberEvents;
}


