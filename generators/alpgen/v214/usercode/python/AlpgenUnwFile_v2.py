import os,sys,logging,gzip
logger = logging.getLogger(__name__)
import AlpgenParameters

class AlpgenUnwFile:
   def __init__(self,filename = None,bytes_per_event = 0):
      self.filename = filename
      self.nevts = 0
      self.bytes_per_event = bytes_per_event
      self.file = None
      self.events_read = 0
      self.current_event = ''
      if filename is not None:
         self.open()

   def get_event(self,eventNumber=None): # event number starting with 0
      
      if self.file is None or self.file.closed:
         logger.error('File is not opened. Get Event failed.')
         return None

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

   def open(self):
      if self.file is not None:
         logger.error(' file already open: ' + self.file.name)
         return -1
      if os.path.exists(self.filename):
         #logger.debug('File opened for reading: ' + self.filename)
         statinfo = os.stat(self.filename)
         self.nevts = statinfo.st_size/self.bytes_per_event
         if '.gz' in self.filename:
            self.file = gzip.open(self.filename,'r')
         else:
            self.file = open(self.filename,'r')
      else:
         #logger.debug('File opened for writing: ' + self.filename)
         self.nevts = 0
         if '.gz' in self.filename:
            self.file = gzip.open(self.filename,'w')
         else:
            self.file = open(self.filename,'w')

      self.events_read = 0
      self.current_event = ''

   def write(self,data):

      if self.file is None or self.file.closed:
         logger.error('File is not opened. Get Event failed.')
         raise Exception('File is not opened. Get Event failed.')

      try:
         self.file.write(data)
      except:
         logger.exception('Error writing to file: ' + str(self.filename))
         raise






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

   process = parfile.get_process()

   unwfile = AlpgenUnwFile(unw_filename,)

   for i in range(options.total_events):
      evt = unwfile.get_event(options.start_event + i)
      if len(evt) > 0:
         print evt,




