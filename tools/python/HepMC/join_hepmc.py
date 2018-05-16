#!/usr/bin/env python
import os,sys,subprocess,optparse,logging,glob,math
logging.basicConfig(level=logging.INFO)
logger=logging.getLogger(__name__)

EVENT_START = 'E '
HEPMC_VERSION_INFO = 'HepMC::Version '
START_OF_EVENT_RECORD = 'HepMC::IO_GenEvent-START_EVENT_LISTING'
END_OF_EVENT_RECORD = 'HepMC::IO_GenEvent-END_EVENT_LISTING'

OUTPUT_FILE_POSTFIX = '%08i.hepmc'

STATUS_INTERVAL=500 # print per events processed

def main():
   parser = optparse.OptionParser(description='join HepMC files into larger files')
   parser.add_option('-g','--input-glob',dest='input_glob',help='glob string of files to join. For example "path/to/files/*.hepmc". Must use the quotations on the command line.')
   parser.add_option('-o','--output-base',dest='output_base',help='The base of the output files. For example, "output_" would produce files "output_0000.hepmc". ',default='output_')
   parser.add_option('-n','--evt-per-file',dest='evt_per_file',help='The number of events to include in each file',type='int',default=5000)
   options,args = parser.parse_args()
   
   if options.input_glob is None:
      parser.error('Must specify -g')

   input_filelist =glob.glob(options.input_glob)
   logger.info(' using glob: ' + options.input_glob + ' found ' + str(len(input_filelist)) + ' files ')
   input_filecount = 0

   output_filename = options.output_base + OUTPUT_FILE_POSTFIX
   logger.info(' using filename pattern: ' + output_filename )
   logger.info(' using ' + str(options.evt_per_file) + ' events per output file.')
   output_filecount = 0
   output_eventcount = 0

   total_eventcount = 0
   
      # extract the header information
   HepMC_Version,StartEvtListing = get_hepmc_header(input_filelist[input_filecount])
   
   # write header to new file
   output_file = open_empty_hepmc_file(output_filename,output_filecount,HepMC_Version,StartEvtListing)

   # loop over all input files and copy to new files
   for input_filename in input_filelist:
      check_hepmc_header(input_filename,HepMC_Version,StartEvtListing)
      input_file = open(input_filename)
      skip_hepmc_header(input_file)
      
      # loop over the lines in the input file
      for line in input_file:
         # reached end of event record
         if line.find(END_OF_EVENT_RECORD) >= 0:
            break
         elif line.find(EVENT_START) >= 0:
            # if output file is full, close it and open a new one
            if output_eventcount >= options.evt_per_file:
               write_hepmc_trailer(output_file)
               output_file.close()
               output_filecount += 1
               output_file = open_empty_hepmc_file(output_filename,output_filecount,HepMC_Version,StartEvtListing)
               output_eventcount = 0
            
            # Event header contains the event number starting at 0 in each file
            # so this needs to be kept straight in the output files
            line = update_event_header(line,output_eventcount)
            
            # logger.info(' processed ' + str(total_eventcount) + ' so far.')
            # increment counters
            output_eventcount += 1
            total_eventcount  += 1
         
         output_file.write(line)
      
      input_file.close()
   # finish off the last file
   write_hepmc_trailer(output_file)
   output_file.close()
   output_filecount += 1

   logger.info(' wrote ' + str(total_eventcount) + ' events to ' + str(output_filecount) + ' files.')


   

# write header information
def write_hepmc_header(file,HepMC_Version,StartEventListing):
   file.write('\n') # first line is always a space
   file.write(HepMC_Version)
   file.write(StartEventListing)

# write trailer information
def write_hepmc_trailer(file):
   file.write(END_OF_EVENT_RECORD + '\n')

def get_output_filename(output_base,filelist):
   nfiles = len(filelist)
   decimals = math.floor(math.log10(nfiles)) + 1
   return output_base + '%0' + str(decimals) + 'i.hepmc'

def get_hepmc_header(filename):
   # open first file
   input_file = open(filename)
   line = ''
   # extract header lines
   # first line is empty
   input_file.readline()
   # second line is version information
   HepMC_Version = input_file.readline()
   if HepMC_Version.find(HEPMC_VERSION_INFO) == -1:
      print 'INVALID FORMAT: Did not find HepMC::Vesion information on correct line'
      print '   instead found: ', HepMC_Version
      return None,None
   # third line is the "Start Event Listing"
   StartEvtListing = input_file.readline()
   if StartEvtListing.find(START_OF_EVENT_RECORD) == -1:
      print 'INVALID FORMAT: Did not find "START_EVENT_LISTING" on correct line'
      print '     instead found: ', StartEvtListing
      return None,None
   
   input_file.close()
   return HepMC_Version,StartEvtListing

def check_hepmc_header(filename,HepMC_Version,StartEvtListing):
   new_vers,new_listing = get_hepmc_header(filename)
   if new_vers != HepMC_Version:
      print 'HEPMC Version difference found.'
   if new_listing != StartEvtListing:
      print 'Event Listing header difference found.'

def skip_hepmc_header(file):
   file.readline() # newline
   file.readline() # Version info
   file.readline() # Event Header info

def open_empty_hepmc_file(filename,count,HepMC_Version,StartEvtListing):
   file = open((filename % count),'w')
   write_hepmc_header(file,HepMC_Version,StartEvtListing)
   return file

def update_event_header(event_header,event_count):
   # split line
   split_line = event_header[0:-1].split()
   # create new event line:
   new_line = split_line[0]
   new_line += ' ' + str(event_count)
   for i in range(2,len(split_line)):
      new_line += ' ' + split_line[i]
   return new_line


if __name__ == '__main__':
   sys.exit(main())

