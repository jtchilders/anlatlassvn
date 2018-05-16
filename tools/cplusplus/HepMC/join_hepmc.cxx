//////////////////////////////////////////////////////////////////////////
// Matt.Dobbs@Cern.CH, Feb 2000
// Example of applying an event selection to the events written to file
// using example_MyPythia.cxx
// Events containing a photon of pT > 25 GeV pass the selection and are
// written to "example_EventSelection.dat"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"
#include "tclap/ArgException.h"

#include <string>
#include <fstream>
#include <iostream>

int main(int argc,char** argv) 
{
   std::string infilelist_filename,outfilename;
   try{
      TCLAP::CmdLine cmd(" Join multiple HepMC files into one ",' ',"0.1");
      TCLAP::ValueArg<std::string> infilelist("i","infilelist","Input file where each line has the name of an input file to be joined",true,"","string",cmd);
      TCLAP::ValueArg<std::string> outfile("o","outfilename","Output filename",false,"output.hepmc","string",cmd);
      cmd.parse(argc,argv);
      infilelist_filename = infilelist.getValue();
      outfilename = outfile.getValue();
   }
   catch (TCLAP::ArgException &e){
      std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
   }
   
   // create output stream
   HepMC::IO_GenEvent ascii_out(outfilename.c_str(),std::ios::out);
   
   // loop over the file list and add each file to the output
   std::ifstream filelist(infilelist_filename.c_str());
   int file_counter = 0;
   int evnt_counter = 0;
   while(filelist.good()){
      char line[500];
      filelist.getline(line,500);
      // open file for input
      std::cout << "Opening file number " << file_counter << " named: " << line << "\n";
      HepMC::IO_GenEvent ascii_in(line,std::ios::in);

      HepMC::GenEvent* evt = ascii_in.read_next_event();
      while (evt){
         if(evnt_counter % 1000 == 0)
            std::cout << "Events Processed: " << evnt_counter << "\n";
         // output event to file
         ascii_out << evt;
         // delete event
         delete evt;
         // read next event
         ascii_in >> evt;
         evnt_counter++;
      }
      file_counter++;
   }

   std::cout << " Processed " << file_counter << " files and " << evnt_counter << " events.\n";
   return 0;
   

}






