#!/usr/bin/env python
import os,sys,subprocess,optparse

EVENT_START = 'E '
HEPMC_VERSION_INFO = 'HepMC::Version '
START_OF_EVENT_RECORD = 'HepMC::IO_GenEvent-START_EVENT_LISTING'
END_OF_EVENT_RECORD = 'HepMC::IO_GenEvent-END_EVENT_LISTING'

STATUS_INTERVAL=500 # print per events processed

def main():
   parser = optparse.OptionParser(description='split HepMC files into smaller chunchs, just specificy events per file.')
   parser.add_option('-n','--evnt-per-file',dest='nevt',help='Number of events per file',type='int')
   parser.add_option('-i','--input-file',dest='input_filename',help='Input filename')
   parser.add_option('-o','--output-base',dest='output_base',help='Output filename base',default='split')
   options,args = parser.parse_args()
   
   if options.nevt is None:
      parser.error('Must supply number of events per file.')
   elif options.input_filename is None:
      parser.error('Must specify input filename')

   try:
      file = open(options.input_filename)
   except:
      print 'Exception while trying to open',options.input_filename,':',sys.exc_info()[1]
      return -1
   
   # extract header lines
   # first line is empty
   file.readline()
   # second line is version information
   HepMC_Version = file.readline()
   if HepMC_Version.find(HEPMC_VERSION_INFO) == -1:
      print 'INVALID FORMAT: Did not find HepMC::Vesion information on correct line'
      print '   instead found: ', HepMC_Version
      return -1
   # third line is the "Start Event Listing"
   StartEvtListing = file.readline()
   if StartEvtListing.find(START_OF_EVENT_RECORD) == -1:
      print 'INVALID FORMAT: Did not find "START_EVENT_LISTING" on correct line'
      print '     instead found: ', StartEvtListing
      return -1


   event_counter = 0
   out_file = None
   file_counter = 0
   for line in file:
      if line.find(EVENT_START) >= 0: # start of new event
         event_counter += 1 # increment event counter
         
         # print status
         if event_counter % STATUS_INTERVAL == 0:
            print 'events processed: ',(file_counter-1)*options.nevt + event_counter,' files_output: ', file_counter

         # close output file if we've reached the max number of events for this file
         if event_counter > options.nevt:
            write_hepmc_trailer(out_file)
            out_file.close()
            out_file = None
            event_counter = 1

         # if no file is open, then open one
         if out_file is None:
            out_file = open(options.output_base + '.' + ('%08i' % file_counter) + '.hepmc','w')
            file_counter += 1 # increment file_counter
            write_hepmc_header(out_file,HepMC_Version,StartEvtListing)
      elif line.find(END_OF_EVENT_RECORD) >= 0: # end of input file
         write_hepmc_trailer(out_file)
         out_file.close()
         out_file = None
         break
      
      out_file.write(line)

   file.close()

   return 0


# write header information
def write_hepmc_header(file,HepMC_Version,StartEventListing):
   file.write('\n') # first line is always a space
   file.write(HepMC_Version)
   file.write(StartEventListing)

# write trailer information
def write_hepmc_trailer(file):
   file.write(END_OF_EVENT_RECORD + '\n')

      



if __name__ == '__main__':
   sys.exit(main())

