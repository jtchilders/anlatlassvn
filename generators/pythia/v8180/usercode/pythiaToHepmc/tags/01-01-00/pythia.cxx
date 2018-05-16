

#include "Pythia8/Pythia.h"
#include "CombineMatchingInput.h"
#include "Pythia8/Pythia8ToHepMC.h"
#include "Gzip_Stream.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"
#include "tclap/ArgException.h"


#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include "sys/stat.h"

#include "mpi.h"

using namespace Pythia8;

#define FINISHED_MSG ">>FINISHED<<\0"
#define HEPMC_EVENT_SIZE_MAX 1000000

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

void dump_event(const unsigned int,const unsigned int,const char*);

#define GZIP_COMPRESSION_LEVEL 1

int main(int argc,char** argv){

   
   // Initialize the MPI environment
   MPI_Init(NULL, NULL);
   // Get the number of processes
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   // Get the rank of the process
   int world_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
   std::stringstream ss_worldrank;
   ss_worldrank.width(8);
   ss_worldrank.fill('0');
   ss_worldrank << world_rank;
   std::string str_world_rank = ss_worldrank.str();
   
   // setup error handling
   MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

   // Get the name of the processor
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   int name_len;
   MPI_Get_processor_name(processor_name, &name_len);
   // Print off a hello world message
   std::cout << "Pythia8 MPI from processor " << processor_name << ", rank " << world_rank
             << " out of " << world_size << "  processors\n";
  
   // parse command line
   std::string outputHepMc,xmldoc,inputPythiaConfig;
   std::vector<std::string> processes;
   unsigned int randomSeed = 1;
   unsigned int numEvts = 0;
   unsigned int ranksPerCollector = 32;
   bool redirectOutput = false;
   bool gzipHepMc = false;
   bool disableOutput = false;
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

      TCLAP::ValueArg<std::string> inputPythiaConfigArg("c","pythia-config","Configuraiton file for Pythia8",true,"pythia.cmnd","string",cmd);

      TCLAP::ValueArg<unsigned int> ranksPerCollectorArg("a","collector","Number of Ranks per Collector. This determines for every N ranks, there will be 1 collector rank which collected the date for the other N-1 ranks (called worker ranks).",false,32,"unsigned int",cmd);
      
      //TCLAP::MultiArg<std::string> processesArg("p","process","Physical Process to generate given in double quotation, acceptable strings can be found in the Pythia8 Documentation, i.e. \"SoftQCD:all\". More than one process can be enabled with multiple calls to this option.",true,"string",cmd);
      
      TCLAP::ValueArg<unsigned int> numEvtsArg( "n","num-evts","Number of events per rank to generate",true,0,"unsigned int",cmd);
      TCLAP::ValueArg<unsigned int> randomSeedArg( "s","random-seed","Random Number Seed, which will have the rank number added to it for each rank. The user should ensure random seeds are not reused to avoid duplication.",false,1,"unsigned int",cmd);

      TCLAP::SwitchArg redirectSwitch("r","redirect-to-file","Redirect stdout and stderr to files.",cmd,false);
      TCLAP::SwitchArg gzipSwitch("g","gzip-output","run gzip on event in memory before writing to file.",cmd,false);

      TCLAP::SwitchArg noOutputSwitch("d","disable-output","do not write any output files (for testing purposes)",cmd,false);
      
      cmd.parse(argc,argv);

      inputPythiaConfig = inputPythiaConfigArg.getValue();
      xmldoc = xmldocArg.getValue();
      //processes = processesArg.getValue();
      numEvts = numEvtsArg.getValue();
      randomSeed = randomSeedArg.getValue();
      ranksPerCollector = ranksPerCollectorArg.getValue();
      redirectOutput = redirectSwitch.getValue();
      gzipHepMc = gzipSwitch.getValue();
      disableOutput = noOutputSwitch.getValue();
      

      std::stringstream filename;
      filename << outputHepMcArg.getValue() << ".rank" << str_world_rank << ".hepmc";
      if(gzipHepMc) filename << ".gz";
      outputHepMc = filename.str();
      
      // check that xmldoc path is correct
      struct stat sb;
      if( stat(xmldoc.c_str(),&sb) != 0 ||  !S_ISDIR(sb.st_mode)){
         std::cerr << " XML config path for Pythia is not set properly. Currently pointing to: " << xmldoc << "\n";
         return -1;
      }

   }
   catch(TCLAP::ArgException &e)  // catch any exceptions
   { std::cerr << "Rank " << str_world_rank << " error: " << e.error() << " for arg " << e.argId() << std::endl; return -1;}

 
   // currently cannot run without at least two ranks so make this check
   if( world_size <= 1 ){
      std::cerr << " Running with only one rank is current not supported, exiting.\n";
      return -1;
   }



   // is this a worker rank or a collector rank?
   bool worker_rank = false, collector_rank = false;
   if(world_rank % ranksPerCollector == 0){
      std::cout << "Rank " << str_world_rank << " is a collector rank.\n";
      collector_rank = true;
   }
   else{
      std::cout << "Rank " << str_world_rank << " is a worker rank.\n";
      worker_rank = true;
   }

   // if this collector rank has no worker ranks, then exit
   if(collector_rank){
      if( world_rank == (world_size - 1) ){
         std::cout << " Rank " << str_world_rank << " has no worker ranks and will exit now.\n";
         MPI_Finalize();
         return 0;
      }
   }




   std::fstream cout_file,cerr_file;
   if(redirectOutput){
      ///////////////////////////////////////////////
      // redirect std::cout std::cerr to files
      std::streambuf *file_cout_buf, *file_cerr_buf, *cout_buf, *cerr_buf;
      // backup cout/cerr stream buffers
      cout_buf = std::cout.rdbuf();
      cerr_buf = std::cerr.rdbuf();
      // open output files
      std::string filename_cout = "stdout.rank" + str_world_rank + ".txt";
      std::string filename_cerr = "stderr.rank" + str_world_rank + ".txt";
      cout_file.open(filename_cout.c_str(),std::fstream::out);
      cerr_file.open(filename_cerr.c_str(),std::fstream::out);
      // get files stream buffer
      file_cout_buf = cout_file.rdbuf();
      file_cerr_buf = cerr_file.rdbuf();
      // assign the file buffers to cout/cerr
      std::cout.rdbuf(file_cout_buf);
      std::cerr.rdbuf(file_cerr_buf);
      /////////////////////////////////////////
   }


   std::cout << "======================================================\n";
   std::cout << " Running Pythia:\n";
   std::cout << "======================================================\n";
   if(collector_rank){
      std::cout << " output filename: " << outputHepMc << "\n";
      std::cout << "======================================================\n";
   }
   else if(worker_rank){
      std::cout << " number of events: " << numEvts << "\n";
      std::cout << "======================================================\n";
      std::cout << " pythia config input from file: " << inputPythiaConfig << "\n";
      std::ifstream file(inputPythiaConfig.c_str());
      char line[500];
      while(!file.eof()){
         file.getline(line,500);
         if(file.eof()) continue;
         std::cout << line << "\n";
      }
      file.close();
      std::cout << "======================================================\n";
   }
   
   // Interface for conversion from Pythia8::Event to HepMC event.
   HepMC::Pythia8ToHepMC* ToHepMC = 0;
   if(worker_rank)
      ToHepMC = new HepMC::Pythia8ToHepMC();

   // Specify file where HepMC events will be stored.
   HepMC::IO_GenEvent* hepMcOutput = 0;

   // keeping this around as example of making gzip stream
   //ogzstream* gzip_outputfile = 0;
   //if(gzipHepMc){
   //   gzip_outputfile = new ogzstream(outputHepMc.c_str(), std::ios::out, ZGIP_COMPRESSION_LEVEL);
   //   hepMcOutput = new HepMC::IO_GenEvent(*gzip_outputfile);
   //}
   //else{
   //   if( world_rank != 0)
   //      hepMcOutput = new HepMC::IO_GenEvent(outputHepMc);
   //}
   
   // This is used by the worker ranks to stream the HepMC event to a stringstream instead of a file stream
   std::stringstream eventString;
   HepMC::IO_GenEvent* sHepMcOutput = 0;
   if(worker_rank)
      sHepMcOutput = new HepMC::IO_GenEvent((std::ostream&)eventString);

   unsigned int nEvent = numEvts;
   unsigned int nAbort = 0;
   
   // Generator and read in commands.
   Pythia* pythia = 0;
   if(worker_rank){
      // create pythia object
      pythia = new Pythia(xmldoc.c_str());
      // load all config options from file
      pythia->readFile(inputPythiaConfig);
   
      // Update Pythia configuration from command line args

      std::stringstream ss;
      ss << "Main:numberOfEvents = " << numEvts;
      pythia->readString(ss.str());
      pythia->readString("Main:timesAllowErrors = 3");
      //pythia->readString("Main:spareMode1 = 0");
      pythia->readString("Init:showChangedSettings = on");
      pythia->readString("Init:showChangedParticleData = on");
      //pythia->readString("Init:showAllSettings = on");
      pythia->readString("Next:numberCount = 500");

      // set the random seed
      {
         unsigned int seed = randomSeed + world_rank;
         std::cout << " Rank " << str_world_rank << " has random seed = " << seed << "\n";
         std::stringstream s_randomSeed; s_randomSeed << "Random:seed = " << seed;
         pythia->readString("Random:setSeed = on");
         pythia->readString(s_randomSeed.str());
      }

      // Initialise Pythia.
      if ( !pythia->init()) {
         std::cerr << "Rank " << str_world_rank << " Error: could not initialise Pythia" << endl;
         return 1;
      }

      // Extract settings to be used in the main program.
      nAbort = pythia->mode("Main:timesAllowErrors");
   }

   // determine message tag
   const unsigned int message_tag = (unsigned int)((float)world_rank/(float)ranksPerCollector) * ranksPerCollector;
   
   ////////////////////////////////////////////////////////////
   // Begin event loop.
   if(worker_rank){ // This is a production rank
      std::cout << " Worker Rank " << str_world_rank << " beginning generation of " << nEvent << " events.\n";
      std::cout << " Worker Rank " << str_world_rank << " sending data to Collector Rank " << message_tag << " \n";

      // this is used for the mpi non-blocking send
      MPI_Request request;
      std::string eventString_copy;

      unsigned int iAbort = 0;
      for (unsigned int iEvent = 0; iEvent < nEvent ;  ++iEvent) {
         //std::cout << " Rank " << str_world_rank << " generating event " << iEvent << "\n";

         // Generate events. Quit if at end of file or many failures.
         if (!pythia->next()) {
            if (++iAbort < nAbort){
               std::cerr << " Worker Rank " << str_world_rank << " event had errors, skipping.\n";
               iEvent--;
               continue;
            }
            std::cerr << " Worker Rank " << str_world_rank << " Abort: too many errors in generation" << endl;
            break;
         }

         // Construct new empty HepMC event and fill it.
         HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
         
         // fill HepMC event with Pythia generated event
         ToHepMC->fill_next_event( *pythia, hepmcevt );

         // before we write a new event to the eventString string stream, we need to make sure the previous event
         // has been sent so that the data is not overwritten
         if(iEvent > 0){
            int flag = 0;
            MPI_Status status;
            MPI_Test(&request,&flag,&status);
            //std::cout << "test flag = " << flag << "\n";
            // flag is set to true if copleted correctly
            if(!flag){
               // sending data is not complete so wait until done
               //std::cout << "Rank " << str_world_rank << " is waiting.\n";
               MPI_Wait(&request,&status);
               //std::cout << "Rank " << str_world_rank << " is done waiting.\n"; 
            }
         
            // reset event string so we don't transmit old information
            eventString.str("");
         

         }

         // write event to stringstream
         (*sHepMcOutput) << hepmcevt;
         
         // check event size, if too big, skip it
         eventString_copy = eventString.str();
         const unsigned int event_size = eventString_copy.size();
         if(event_size > HEPMC_EVENT_SIZE_MAX ){
            std::cerr << " Worker Rank " << str_world_rank << " Error: Event size is too large, cannot send over MPI so skipping it. Size = " << event_size << "\n";
            continue;
         }
         // std::cout << " Rank " << str_world_rank << " Event " << iEvent << " sending " << event_size << " characters \n";
         // transmit event string to head rank
         //dump_event(world_rank,iEvent,eventString_copy.c_str());
         int err = MPI_Isend(&eventString_copy[0],event_size,MPI_CHAR,message_tag,message_tag,MPI_COMM_WORLD,&request);
         if( err != MPI_SUCCESS){
            std::cerr << " Worker Rank " << str_world_rank << " Error sending event.\n";
         }
         

         delete hepmcevt;
         hepmcevt = 0;
         
         // store events until max_store_event_count have been collected
         //if(event_store.size() < max_store_event_count){
         //   event_store.push_back(hepmcevt);
         //}
         //else{ // then write them all at once
            //for(unsigned int store_i=0;store_i<event_store.size();++store_i){
               // Write the HepMC event to file. Done with it.
               //if(!disableOutput)
                  //(*hepMcOutput) << hepmcevt; //event_store[store_i];
               //delete hepmcevt; //event_store[store_i];
               //event_store[store_i] = 0;
            //}// end for
            // clear event store to reset
            //event_store.clear();
         //} // end else
         // End of event loop.
      }
      // send message to lead rank that finished event generation
      
      int flag = 0;
      MPI_Status status;
      MPI_Test(&request,&flag,&status);
      // flag is set to true if copleted correctly
      if(!flag){
         // sending data is not complete so wait until done
         //std::cout << "Rank " << str_world_rank << " is waiting.\n";
         MPI_Wait(&request,&status);
      }
    
      eventString.str("");
      eventString << FINISHED_MSG;
      MPI_Send(&eventString.str()[0],eventString.str().size(),MPI_CHAR,message_tag,message_tag,MPI_COMM_WORLD);
      std::cout << " Worker Rank " << str_world_rank << " has finished generating events.\n";
   }
   else{ // rank 0 retreives all the data and writes it to disk
      MPI_Status status;
      // create buffer to receive event
      const unsigned int buffer_size = HEPMC_EVENT_SIZE_MAX;
      char buffer[buffer_size];
      unsigned int finished_messages_received = 0;
      unsigned int finished_messages_expected = world_size - 1;
      // if the world size is bigger than the number of ranks per collector
      // resize the expected number of ranks to get messages from
      if(world_size > ranksPerCollector){
         finished_messages_expected = ranksPerCollector - 1;
      }
      std::cout << "Collector Rank " << str_world_rank << " collecting data from " << finished_messages_expected << " production ranks.\n";
      // open outputfile
      //ogzstream* gzip_outputfile = 0;
      std::ostream* output = 0;
      if(gzipHepMc){
         output = new ogzstream(outputHepMc.c_str(), std::ios::out, GZIP_COMPRESSION_LEVEL);
      }
      else{
         output = new ofstream(outputHepMc.c_str());
      }

      unsigned int nevts = 0;
      // loop over the expected number of events
      while(true){
         //MPI_Request request;
         int err = MPI_Recv(buffer,buffer_size,MPI_CHAR,MPI_ANY_SOURCE,message_tag,MPI_COMM_WORLD,&status);
         if(err != MPI_SUCCESS){
            std::cerr << "Collector Rank " << str_world_rank << " Error receiving event.\n";
            continue;
         }
         //int err = MPI_Wait(&request,&status);
         //std::cout << " status = " << err << " " << MPI_SUCCESS << "\n";
         // get how many characters were received
         int recv_count = 0;
         MPI_Get_count(&status,MPI_CHAR,&recv_count);
         //std::cout << " Rank " << str_world_rank << " Event " << nevts << " received " << recv_count << " characters\n";
         if( (unsigned int)recv_count >= buffer_size){
            std::cerr << "Collector Rank " << str_world_rank << " has a character count of " << recv_count << " which is larger than the allowed buffer size of " << buffer_size << ", skipping event.\n";
            continue;
         }
         // must null terminate the string so that old data is not dumped to file.
         //std::cout << " recv_count = " << recv_count << "\n";
         buffer[recv_count] = '\0';
         //std::cout << buffer << "\n\n";
         //dump_event(world_rank,nevts,buffer);
         if( strcmp(buffer,FINISHED_MSG) == 0){
            std::cout << "Collector Rank " << str_world_rank << " received FINISHED_MSG from rank " << status.MPI_SOURCE << "\n";
            //std::cout << " buffer = " << buffer << "\n";
            finished_messages_received++;
            if( finished_messages_received >= finished_messages_expected ){
               std::cout << "Collector Rank " << str_world_rank << " finished collecting events.\n";
               break;
            }
         }
         (*output) << buffer;
         nevts++;

         if( nevts % 500 == 0)
            std::cout << "Collector Rank " << str_world_rank << " received " << nevts << " events\n";
      }
      // Add HepMC trailer
      (*output) << "HepMC::IO_GenEvent-END_EVENT_LISTING";
      //output->close();

      delete output;
      output =0;
   }
   // Final statistics and done.
   if(pythia) pythia->stat();

   if(redirectOutput){
      cout_file.close();
      cerr_file.close();
   }

   //if(gzipHepMc){
   //   gzip_outputfile->close();
   //   delete gzip_outputfile;
   //   gzip_outputfile = 0;
   //}
   
   if( hepMcOutput){
      delete hepMcOutput;
      hepMcOutput = 0;
   }
   if( sHepMcOutput){
      delete sHepMcOutput;
      sHepMcOutput = 0;
   }

   if(ToHepMC){
      delete ToHepMC;
      ToHepMC = 0;
   }

   // Finalize the MPI environment.
   MPI_Finalize();

   std::cout << " Rank " << str_world_rank << " exiting.\n";



   return 0;
}

void dump_event(const unsigned int rank,const unsigned int evtnum,const char* buffer){
   std::stringstream ss;
   ss << "rank" << rank << "_event" << evtnum << ".txt";
   ofstream out(ss.str().c_str());
   out << buffer;
   out.close();
}


