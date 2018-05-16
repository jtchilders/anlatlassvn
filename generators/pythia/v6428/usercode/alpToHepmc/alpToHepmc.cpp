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

#include "sys/stat.h"

#include "HepMC/PythiaWrapper.h"
#include "HepMC/IO_HEPEVT.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/IO_AsciiParticles.h"
#include "HepMC/GenEvent.h"
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"
#include "tclap/ArgException.h"

/*
c general parameters
      integer nparam
      parameter (nparam=600)
      integer parlen,partyp
      character chpar*8,chpdes*70
      double precision parval
      common/AHpars/parval(nparam),chpar(nparam),chpdes(nparam)
     $     ,parlen(nparam),partyp(nparam)
*/

const int npars = 600;
extern "C" {
   extern struct{
      double parval[npars];
      char chpar[npars];
      char chpdes[npars];
      int parlen[npars];
      int partyp[npars];
   }ahpars_;
}
#define ahpars_vars ahpars_;
/*
      CHARACTER CBASENAME*256
      COMMON/ALPGEN/CBASENAME
*/
const int FORTRAN_STRING_SIZE=256;
extern "C" {
   extern struct {
      char cbasename[FORTRAN_STRING_SIZE];
   } alpgen_;
}
#define alpgen_vars alpgen_

int get_num_events(std::string filename);

int main(int argc,char** argv){

   // parse command line
   std::string outputHepMc,inputAlpgenBase;
   int maxevts = -1;
   bool jetInclusive = false;
   try{
      TCLAP::CmdLine cmd("Hadronize Alpgen events using Pythia8. Takes Alpgen *.unw and *_unw.par files and outputs HepMC format.", ' ', "0.1");
      
      
      TCLAP::ValueArg<std::string> outputHepMcArg("o","output-file","Output Filename",false,"pythia_output.hepmc","string",cmd);

       TCLAP::ValueArg<std::string> inputAlpgenBaseArg("a","alpgen-file","Input Filename base for Alpgen inputs. 'alpout' for 'alpout.unw' and 'alpout_unw.par'.",true,"alpout","string",cmd);

      TCLAP::ValueArg<int> maxevtsArg("n","numevts","Number of events to process. If not specified or number is negative, all events will be processed",true,-1,"int",cmd);
      TCLAP::SwitchArg jetInclArg("i","inclusive","Jet matching is by default exclusive, but this flag makes it inclusive",cmd,false);
      
      cmd.parse(argc,argv);

      outputHepMc = outputHepMcArg.getValue();
      inputAlpgenBase = inputAlpgenBaseArg.getValue();
      maxevts = maxevtsArg.getValue();
      jetInclusive = jetInclArg.getValue();
      

      std::cout << " Running with settings: \n";
      std::cout << "  input AlpGen base:  " << inputAlpgenBase << "\n";
      std::cout << "  ouput HepMC file:   " << outputHepMc << "\n";
      if(jetInclusive)
         std::cout << "  using inclusive jet matching.\n";
      else
         std::cout << "  using exclusive jet matching.\n";
      if(maxevts < 0)
         std::cout << "  running over all input events\n";
      else
         std::cout << "  running over " << maxevts << " input events\n";

   }
   catch(TCLAP::ArgException &e)  // catch any exceptions
   { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; return -1;}


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

   // set Alpgen Variables
   for(unsigned int i=0;i<FORTRAN_STRING_SIZE;++i)
      if (i<inputAlpgenBase.size())
         alpgen_vars.cbasename[i] = inputAlpgenBase[i];
      else
         alpgen_vars.cbasename[i] = ' ';
   
   // keep in mind that index 500 is 501 in fortran counting
   ahpars_.parval[500] = 20.1;
   ahpars_.parval[501] = 0.41;
   ahpars_.parval[502] = 6.01;
   if(jetInclusive)
      ahpars_.parval[503] = 0;
   else
      ahpars_.parval[503] = 1;

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
   HepMC::IO_GenEvent ascii_io(outputHepMc,std::ios::out);
   
   int eventMax = maxevts;
   if(eventMax < 0){
      std::cout << "ERROR getting number of events for this run\n";
      return -1;
   }

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


