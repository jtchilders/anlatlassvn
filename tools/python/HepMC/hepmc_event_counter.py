#!/usr/bin/env python
import os,sys,optparse,glob

EVENT_START = 'E '

def main():
   parser = optparse.OptionParser(description='Count the number of events in a HepMC file.')
   parser.add_option('-i','--input-file',dest='input_filename',help='Input filename')
   parser.add_option('-g','--input-glob',dest='input_glob',help='Input glob string for looping over many files, must use quotation marks to enclose string pattern.')
   options,args = parser.parse_args()
   
   if options.input_filename is None and options.input_glob is None:
      parser.error('Must specify input filename or input glob')

   event_count = 0

   if options.input_filename is not None:
      event_count = count_events_in_hepmc_file(options.input_filename)

   if options.input_glob is not None:
      event_count += count_events_in_hepmc_files(options.input_glob) 


def count_events_in_hepmc_files(glob_string):
   files = glob.glob(glob_string)
   event_count = 0
   for file in files:
      event_count += count_events_in_hepmc_file(file)

   print ' Total Events = ' + str(event_count)

   return event_count

def count_events_in_hepmc_file(filename):
   try:
      file = open(filename)
   except:
      print 'Exception while trying to open',filename,':',sys.exc_info()[1]
      return -1
  
   
   event_count = 0
   for line in file:
      if line.find(EVENT_START) >= 0:
         event_count += 1

   file.close()
   
   print 'Number of Events in file',filename,' = ',event_count
   
   return event_count



if __name__ == '__main__':
   main()

