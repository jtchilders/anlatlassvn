

#include "Pythia8/Pythia.h"
#include "CombineMatchingInput.h"
#include "Pythia8/Pythia8ToHepMC.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"
#include "tclap/ArgException.h"


#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include "sys/stat.h"

#include "mpi.h"

using namespace Pythia8;

class PythiaProcess{
public:
   PythiaProcess(){};
   ~PythiaProcess(){};
   
   enum Processes{
      SOFTQCD_ALL=0,
      SOFTQCD_NONDIFFRACTIVE,
      SOFTQCD_ELASTIC,
      SOFTQCD_SINGLEDIFFRACTIVE,
      SOFTQCD_DOUBLEDIFFRACTIVE,
      SOFTQCD_CENTRALDIFFRACTIVE,
      SOFTQCD_INELASTIC,
      
      HARDQCD_ALL,
      HARDQCD_GG2GG,
      HARDQCD_GG2QQBAR,
      HARDQCD_QG2QG,
      HARDQCD_QQ2QQ,
      HARDQCD_QQBAR2GG,
      HARDQCD_QQBAR2QQBARNEW,
      HARDQCD_NQUARKNEW,

      HARDQCD_GG2CCBAR,
      HARDQCD_QQBAR2CCBAR
   };


};


int main(int argc,char** argv){
  
   // Initialize the MPI environment
   MPI_Init(NULL, NULL);
   // Get the number of processes
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   // Get the rank of the process
   int world_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

   // Get the name of the processor
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   int name_len;
   MPI_Get_processor_name(processor_name, &name_len);
   // Print off a hello world message
   std::cout << "Pythia8 MPI from processor " << processor_name << ", rank " << world_rank
             << " out of " << world_size << "  processors\n";

   // parse command line
   std::string outputHepMc,xmldoc;
   std::vector<std::string> processes;
   unsigned int numEvts;
   try{
      TCLAP::CmdLine cmd("Run Pythia8.Outputs HepMC format.", ' ', "0.1");
      
      // Define a value argument and add it to the command line.
      // A value arg defines a flag and a type of value that it expects,
      // such as "-o alpgen.hepmc" or "--output-file=alpgen.hepmc".
      // true/false defines if the value is required
      // "pythia_output.hepmc" is a default to use if it is not defined on the command line
      // "string" is the human readable version of the data type expected
      // "cmd" add this argument to the command line object
      TCLAP::ValueArg<std::string> outputHepMcArg("o","output-filebase","Output File base-name, all files will be prefixed with '.hepmc'",false,"pythia_output","string",cmd);

      TCLAP::ValueArg<std::string> xmldocArg("x","xmldoc","Path to the xmldoc/Index.xml file for Pythia setup. Typically this is the installation path of Pythia.",false,"./xmldoc","string",cmd);

      //TCLAP::ValueArg<std::string> inputPythiaConfigArg("c","pythia-config","Configuraiton file for Pythia8",false,"pythia.cmnd","string",cmd);
      
      TCLAP::MultiArg<std::string> processesArg("p","process","Physical Process to generate given in double quotation, acceptable strings can be found in the Pythia8 Documentation, i.e. \"SoftQCD:all\". More than one process can be enabled with multiple calls to this option.",true,"string",cmd);
      
      TCLAP::ValueArg<unsigned int> numEvtsArg( "n","num-evts","Number of events to generate",true,0,"unsigned int",cmd);
      
      cmd.parse(argc,argv);

      //inputPythiaConfig = inputPythiaConfigArg.getValue();
      std::stringstream filename;
      filename << outputHepMcArg.getValue() << ".rank";
      filename.width(5);
      filename.fill('0');
      filename << world_rank;
      filename << ".hepmc";
      outputHepMc = filename.str();
      xmldoc = xmldocArg.getValue();
      processes = processesArg.getValue();
      numEvts = numEvtsArg.getValue();
      // check that xmldoc path is correct
      struct stat sb;
      if( stat(xmldoc.c_str(),&sb) != 0 ||  !S_ISDIR(sb.st_mode)){
         std::cerr << " XML config path for Pythia is not set properly. Currently pointing to: " << xmldoc << "\n";
         return -1;
      }

   }
   catch(TCLAP::ArgException &e)  // catch any exceptions
   { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; return -1;}
  
   std::cout << "===========================\n";
   std::cout << "Generating processes:\n";
   for(unsigned int i=0;i<processes.size();++i)
      std::cout << "     " << i << ". " << processes[i] << "\n";
   std::cout << "===========================\n";
   std::cout << " output filename: " << outputHepMc << "\n";
   std::cout << "===========================\n";
   
   // Interface for conversion from Pythia8::Event to HepMC event.
   HepMC::Pythia8ToHepMC ToHepMC;

   // Specify file where HepMC events will be stored.
   HepMC::IO_GenEvent ascii_io(outputHepMc, std::ios::out);
   
   // Generator and read in commands.
   Pythia pythia(xmldoc.c_str());
   //pythia.readFile(inputPythiaConfig);
   
   // Configure Pythia
   std::stringstream ss;
   ss << "Main:numberOfEvents = " << numEvts;
   pythia.readString(ss.str());
   pythia.readString("Main:timesAllowErrors = 3");
   pythia.readString("Main:spareMode1 = 0");
   pythia.readString("Init:showChangedSettings = on");
   pythia.readString("Init:showChangedParticleData = on");
   //pythia.readString("Init:showAllSettings = on");
   pythia.readString("Next:numberCount = 500");
   // add processes
   for(unsigned int i=0;i<processes.size();++i){
      ss.str("");
      ss << processes[i] << " = on";
      pythia.readString(ss.str());
   }

   // set the random seed
   {
      unsigned int randomSeed = numEvts*world_rank;
      std::stringstream s_randomSeed; s_randomSeed << "Random:seed = " << randomSeed;
      pythia.readString("Random:setSeed = on");
      pythia.readString(s_randomSeed.str());
   }

   // Initialise Pythia.
   if ( !pythia.init()) {
      cout << "Error: could not initialise Pythia" << endl;
      return 1;
   };

   // Extract settings to be used in the main program.
   int nEvent = pythia.mode("Main:numberOfEvents");
   int nAbort = pythia.mode("Main:timesAllowErrors");
   int nSkip  = pythia.mode("Main:spareMode1");

   // Begin event loop. Optionally quit it before end of file.
   int iAbort = 0;
   for (int iEvent = 0; ;  ++iEvent) {
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
   
   // Finalize the MPI environment.
   MPI_Finalize();

   return 0;
}


