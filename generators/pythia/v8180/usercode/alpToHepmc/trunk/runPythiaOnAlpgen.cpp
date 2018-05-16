

#include "Pythia8/Pythia.h"
#include "CombineMatchingInput.h"
#include "Pythia8/Pythia8ToHepMC.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"
#include "tclap/ArgException.h"


#include <string>
#include <iostream>
#include "sys/stat.h"

using namespace Pythia8;

int main(int argc,char** argv){
   
   // parse command line
   std::string inputAlpgenBase,outputHepMc,xmldoc;
   int maxevts = -1;
   bool jetInclusive = false;
   try{
      TCLAP::CmdLine cmd("Hadronize Alpgen events using Pythia8. Takes Alpgen *.unw and *_unw.par files and outputs HepMC format.", ' ', "0.1");
      
      // Define a value argument and add it to the command line.
      // A value arg defines a flag and a type of value that it expects,
      // such as "-a alpgen" or "--alpgen-base=alpgen".
      // true/false defines if the value is required
      // "mc.root" is a default to use if it is not defined on the command line
      // "string" is the human readable version of the data type expected
      // "cmd" add this argument to the command line parser
      TCLAP::ValueArg<std::string> inputAlpgenBaseArg("a","alpgen-base","Base name for the Alpgen input, i.e., to read alpgen data from 'alpgen.unw' and 'alpgen_unw.par' files pass 'alpgen'",false,"alpout","string",cmd);
      
      TCLAP::ValueArg<std::string> outputHepMcArg("o","output-file","Output Filename",false,"pythia_output.hepmc","string",cmd);

      TCLAP::ValueArg<std::string> xmldocArg("x","xmldoc","Path to the xmldoc/Index.xml file for Pythia setup.",false,"../xmldoc","string",cmd);

      TCLAP::ValueArg<int> maxevtsArg("n","numevts","Number of events to process. If not specified or number is negative, all events will be processed",false,-1,"int",cmd);
      TCLAP::SwitchArg jetInclArg("i","inclusive","Jet matching is by default exclusive, but this flag makes it inclusive",cmd,false);
      
      cmd.parse(argc,argv);

      inputAlpgenBase = inputAlpgenBaseArg.getValue();
      outputHepMc = outputHepMcArg.getValue();
      xmldoc = xmldocArg.getValue();
      maxevts = maxevtsArg.getValue();
      jetInclusive = jetInclArg.getValue();
      // check that xmldoc path is correct
      struct stat sb;
      if( stat(xmldoc.c_str(),&sb) != 0 ||  !S_ISDIR(sb.st_mode)){
         std::cerr << " XML config path for Pythia is not set properly. Currently pointing to: " << xmldoc << "\n";
         return -1;
      }

      std::cout << " Running with settings: \n";
      std::cout << "  input AlpGen base:  " << inputAlpgenBase << "\n";
      std::cout << "  ouput HepMC file:   " << outputHepMc << "\n";
      std::cout << "  pythia xmldoc path: " << xmldoc << "\n";
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
   
   
   // Interface for conversion from Pythia8::Event to HepMC event.
   HepMC::Pythia8ToHepMC ToHepMC;

   // Specify file where HepMC events will be stored.
   HepMC::IO_GenEvent ascii_io(outputHepMc, std::ios::out);
   
   // Generator and read in commands.
   Pythia pythia(xmldoc.c_str());
   
   // Configure Pythia
   pythia.readString("Main:numberofEvents = -1");        // -1 for all
   pythia.readString("Main:timesAllowErrors = 3");       // how many aborts before run stops
   pythia.readString("Main:spareMode1 = 0");             // skip n events at beginning of file
   pythia.readString("Init:showChangedSettings = on");   // list changed settings
   pythia.readString("Init:showChangedParticleData = on"); // list changed particle data
   pythia.readString("Next:numberCount = 500");          // print message every n events
   pythia.readString("JetMatching:merge = on");          // turn on jet matching
   pythia.readString("JetMatching:scheme = 2");          // use jet matching scheme for Alpgen (1 = MadGraph version)
   if(jetInclusive)
      pythia.readString("JetMatching:exclusive = 0");    // All partons must match jets, but additional jets are allowed.
   else
      pythia.readString("JetMatching:exclusive = 1");    // All partons must match jets, and no additional jets are allowed.
   pythia.readString("Alpgen:setMLM = on");              // Set JetMatching variables eTjetMin,coneRadius,etaJetMax based on Alpgen input
    
   if(inputAlpgenBase.size() > 0){
      std::cout << " Setting Alpgen input file basename to " << inputAlpgenBase << "\n";
      pythia.readString(std::string("Alpgen:file = ") + inputAlpgenBase); 
   }

   // Extract settings to be used in the main program.
   int nEvent = pythia.mode("Main:numberOfEvents");
   int nAbort = pythia.mode("Main:timesAllowErrors");
   int nSkip  = pythia.mode("Main:spareMode1");

   // Create UserHooks pointer. Stop if it failed. Pass pointer to Pythia.
   CombineMatchingInput combined;
   UserHooks* matching = combined.getHook(pythia);
   if (!matching){
      cout << "Error: could not get hook for CombineMatchingInput.\n";
      return -1;
   }
   pythia.setUserHooksPtr(matching);
   
   //pythia.readString("Init:showAllSettings = on");

   // Initialise Pythia.
   if ( !pythia.init()) {
      cout << "Error: could not initialise Pythia" << endl;
      return 1;
   };

   // Optionally skip ahead in LHEF.
   pythia.LHAeventSkip( nSkip );

   // Begin event loop. Optionally quit it before end of file.
   int iAbort = 0;

   //for (int iEvent = 0; iEvent < maxevts || maxevts < -1 ;  ++iEvent) {
   for (int iEvent = 0; iEvent < maxevts || maxevts <= -1 ;  ++iEvent) {
      if (nEvent > 0 && iEvent >= nEvent) break;

      // Generate events. Quit if at end of file or many failures.
      if (!pythia.next()) {
         if (pythia.info.atEndOfFile()) {
            cout << "Info: end of input file reached" << endl;
            break;
         }
         if (++iAbort < nAbort) continue;
         cout << "Abort: too many errors in generation" << endl;
         break;
      }

      // Event analysis goes here.
      
      // Construct new empty HepMC event and fill it.
      // Units will be as chosen for HepMC build; but can be changed
      // by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      ToHepMC.fill_next_event( pythia, hepmcevt );

      // Write the HepMC event to file. Done with it.
      ascii_io << hepmcevt;
      delete hepmcevt;
      // End of event loop.
   }

   // Final statistics and done.
   pythia.stat();
   delete matching;
   
   return 0;
}


