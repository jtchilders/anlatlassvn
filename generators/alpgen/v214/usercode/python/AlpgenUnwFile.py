import os,sys,logging
logger = logging.getLogger(__name__)
from FileTools import GetLineCount
from AlpgenUnwParFile import AlpgenUnwParFile
import AlpgenParameters
ALPGEN_ATHENA_FILEBASE='alpgen.XXXXXX.TXT.v1._'

class AlpgenUnwFile:
   STATUS_INTERVAL = 500000

   def __init__(self,filename,par_file):
      self.filename = filename
      if not os.path.exists(self.filename):
         logger.error('File not found: ' + self.filename)
         return
      self.parfile = par_file
      self.nevts = par_file.event_count
      self.linecount = GetLineCount(self.filename)
      self.lines_per_event = self.linecount / self.nevts
      self.bytes_per_event = None
      try:
         # find njets
         njets = None
         for parameter in self.parfile.parameters:
            if 'njets' in parameter.comment:
               njets = parameter.value
               break
         if njets is not None:
            self.bytes_per_event = AlpgenParameters.event_sizes[self.parfile.get_process()][njets]
      except:
         logger.exception('Error setting up Alpgen Input File: ' + str(self.filename))
         raise
      self.file = None
      self.events_read = 0
      self.current_event = ''

   def get_event(self,eventNumber=None): # event number starting with 0
      # make sure file is open
      if self.file is None or self.file.closed:
         self.file = open(self.filename)
         self.events_read = 0

      # move file pointer to event requested
      if eventNumber is not None and self.parfile.get_process() is not None:
         try:
            self.file.seek(eventNumber*self.bytes_per_event)
         except:
            logger.exception('Error trying to seek through file for event. EventNum = ' + str(eventNumber) + ' process = ' + str(self.parfile.get_process()))
            raise

      # read event
      if self.parfile.get_process() is not None:
         try:
            self.current_event = self.file.read(self.bytes_per_event)
            if len(self.current_event) < self.bytes_per_event:
               logger.error('Requested ' + str(self.bytes_per_event) + ' bytes from file but only received ' + str(len(self.current_event)) + '. May have unexpectedly reached the end of file.')
               return ''
            else:
               self.events_read += 1
         except:
            logger.error('Error trying to read event from file. EventNum = ' + str(self.events_read) + ' process = ' + str(self.parfile.get_process()))

      return self.current_event





   def SplitFile(self,
                  events_per_file,
                  output_base,
                  max_files = -1,
                  file_number_offset = 0,
                  skip_last_file = False,
                  first_event = 0,
                  first_file = 0,
                 ):
      
      try:
         file = open(self.filename)
      except:
         logger.exception('Error opening file ' + self.filename + '. Exception: ' + sys.exec_info()[1])
         raise
      
      logger.info('\n source file: ' + self.filename + '\n will be split into '
                   + str(self.nevts / events_per_file + 1) + ' files with offset: ' + str(file_number_offset)
                 )

      file_event_counter = 0  # count events that have been written to current file
      total_event_counter = 0 # count total events that have been read
      line_counter  = 0 # count the number of lines that have been read
      file_counter  = 0 # count the actual number of files that have been written, used for file name numbering
      out_file = None
      output_unw_filename = ''
      output_unw_par_filename = ''
      base = ''
      output_basenames = []
      for line in file:
         if line_counter % self.lines_per_event == 0: # start of new event
            total_event_counter += 1 # increment total event counter
            
            # if this event number is larger than the first_event
            # then use this event in the output
            if total_event_counter >= first_event:
               file_event_counter += 1 # increment file event counter

               if file_event_counter % AlpgenUnwFile.STATUS_INTERVAL == 0:
                  logger.info('events processed: ' + str((file_counter-1)*events_per_file + file_event_counter) 
                              + ' files output: ' + str(file_counter) + ' line counter: ' + str(line_counter)
                             )

               # close output file if we've reached the max nubmer of events for this file
               if file_event_counter > events_per_file and out_file is not None:
                  # stop if max_files is reached
                  if max_files >= 0 and file_counter >= max_files:
                     break
                  out_file.close()
                  out_file = None
                  # write new par file for this data file
                  self.parfile.event_count = events_per_file
                  self.parfile.lumi = self.parfile.event_count / self.parfile.cross_section
                  self.parfile.write_file(output_unw_par_filename)
                  

               # if no file is open, then open a new one
               if out_file is None:
                  # if skip_last_file & the number of events left to write are less than what each file should contain, then exit this loop
                  if skip_last_file and (self.nevts - total_event_counter ) < events_per_file:
                     break
                  # check that file count is above first_file or else skip this
                  if first_file <= int( (total_event_counter-1) / events_per_file):
                     base = output_base + ('%05i' % (file_counter+file_number_offset))
                     output_unw_filename = base  + '.unw'
                     output_unw_par_filename = base + '_unw.par'
                     out_file = open(output_unw_filename,'w')
                     output_basenames.append(base)
                     file_counter += 1 # increment file counter
                     file_event_counter = 1 # reset file_event_counter
                  
         # write the line to file
         if out_file is not None:
            out_file.write(line)
         
         line_counter += 1
      
      # reached end of file so close the output file and write a new par file
      if out_file is not None:
         out_file.close()
         out_file = None
         self.parfile.event_count = file_event_counter
         self.parfile.lumi = self.parfile.event_count / self.parfile.cross_section
         self.parfile.write_file(output_unw_par_filename)

      return output_basenames

   def MakeDataSetTarball(self,basenames,tarball_base = ALPGEN_ATHENA_FILEBASE,counter_offset = 0):

      if '._' not in tarball_base[-2:]: tarball_base += '._'

      # loop over basenames
      for i in range(len(basenames)):
         basename = basenames[i]
         unw_filename = basename + '.unw'
         par_filename = basename + '_unw.par'
         
         str_counter = ('%05i' % (i+counter_offset))
         new_unw_filename = os.path.basename(tarball_base) + str_counter + '.events'
         new_par_filename = os.path.basename(tarball_base) + str_counter + '.dat'
         tarball_filename = tarball_base + str_counter + '.tar.gz'

         os.system('cp ' + unw_filename + ' ' + new_unw_filename)
         os.system('cp ' + par_filename + ' ' + new_par_filename)

         os.system('tar zcf ' + tarball_filename + ' ' + new_unw_filename + ' ' + new_par_filename)
         os.system('rm -f ' + new_unw_filename + ' ' + new_par_filename)


if __name__ == '__main__':
   import optparse
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='')
   parser.add_option('-i','--input-basename',dest='input_basename',help='Input basename from which to read. Assumed that basename + .unw are the events and basename + _unw.par is the parameter file')
   parser.add_option('-s','--start-event',dest='start_event',help='Starting reading at this event, counting starts with 0.',default=0,type='int')
   parser.add_option('-n','--total-events',dest='total_events',help='Total number of events to print to screen',type='int',default=1)
   options,args = parser.parse_args()

   
   manditory_args = [
                     'input_basename',
                  ]

   for man in manditory_args:
      if not options.__dict__[man]:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)
   
   unw_filename = options.input_basename + '.unw'
   par_filename = options.input_basename + '_unw.par'

   if not os.path.exists(unw_filename):
      logger.error(' no unweighted events file: ' + unw_filename)
      sys.exit(-1)
   if not os.path.exists(par_filename):
      logger.error(' no unweighted parameters file: ' + par_filename)
      sys.exit(-2)
   
   import AlpgenUnwParFile
   parfile = AlpgenUnwParFile.AlpgenUnwParFile.read_file(par_filename)

   unwfile = AlpgenUnwFile(unw_filename,parfile)

   for i in range(options.total_events):
      evt = unwfile.get_event(options.start_event + i)
      if len(evt) > 0:
         print evt,




