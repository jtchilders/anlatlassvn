import os,sys,logging,shutil
logger = logging.getLogger(__name__)

import AlpgenParameters

class AlpgenUnwParFile:
   EVENT_COUNT_ID='unwtd events, lum (pb-1)'
   CROSS_SECTION_ID='Crosssection +- error'
   EVENT_COUNT_COLUMN=12
   LUMI_COLUMN=33
   BEGIN_PARAMETERS='************** run parameters'
   END_PARAMETERS='************** end parameters'
   COMMENT_CHAR='!'
   HARD_PROCESS_ID='hard process code'
   MASSES_ID='mc,mb,mt,mw,mz,mh'
   EXCLUSIVE_JET_MATCHING=1
   INCLUSIVE_JET_MATCHING=0
   
   

   class AlpgenUnwParameter:
      def __init__(self,key,value,comment):
         self.key       = key
         self.value     = value
         self.comment   = comment
      def get_line(self):
         return '%4i %15.5f ! %s' % (self.key,self.value,self.comment)
      def __str__(self):
         return self.get_line()
      @staticmethod
      def read_line(line):
         split_index = line.find(AlpgenUnwParFile.COMMENT_CHAR)
         split_values = line[0:split_index].split()
         key = int(split_values[0])
         value = float(split_values[1])
         comment = line[split_index+1:]
         return AlpgenUnwParFile.AlpgenUnwParameter(key,value,comment)

   def __init__(self):
      self.hard_process_code = 0
      self.mc = 0
      self.mb = 0
      self.mt = 0
      self.mw = 0
      self.mz = 0
      self.mh = 0
      self.cross_section = 0
      self.cross_section_error = 0
      self.parameters = []
      self.lumi = 0
      self.event_count = 0
      self.filename = None
      self.start_lines = []

   def write_file(self,new_filename):
      if self.filename == new_filename:
         logger.error('Cannot overwrite existing file.')
         return
      try:
         new_file = open(new_filename,'w')
      except:
         logger.exception('Error opening file ' + new_filename + '. Exception: ' + sys.exec_info()[1])
         return
      
      for line in self.start_lines:
         new_file.write(line)

      new_file.write(AlpgenUnwParFile.BEGIN_PARAMETERS + '\n')
      new_file.write('%4i ! hard process code\n' % self.hard_process_code)
      new_file.write(' %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f ! mc,mb,mt,mw,mz,mh\n' % (self.mc,self.mb,self.mt,self.mw,self.mz,self.mh))
      for par in self.parameters:
         new_file.write(par.get_line() + '\n')
      new_file.write(AlpgenUnwParFile.END_PARAMETERS + '\n')

      new_file.write('%25.8f %24.8f  ! Crosssection +- error (pb)\n' % (self.cross_section,self.cross_section_error))
      new_file.write('%10i %18.8f  ! unwtd events, lum (pb-1)\n' % (self.event_count,self.lumi))
      new_file.close()

   def get_process(self):
      if self.hard_process_code != 0:
         try:
            return AlpgenParameters.hard_process_codes[self.hard_process_code]
         except:
            logger.error('Hard Process Code ' + str(self.hard_process_code) + ' not found in dictionary.')
            return None
      else:
         logger.error('Hard Process Code not yet set. Perhaps no file has been read yet?')
      return None

   def get_njets(self):
      njets_key = 10
      njets = self.get_parameter_by_key(njets_key)
      if njets is not None:
         return int(njets)
      return None

   def get_drjmin(self):
      drjmin_key = 50
      drjmin = self.get_parameter_by_key(drjmin_key)
      if drjmin is not None:
         return float('%0.2f' % float(drjmin))
      return None

   def get_parameter_by_key(self,key):
      for parameter in self.parameters:
         if parameter.key == key:
            return parameter.value
      return None

   def get_bytes_per_event(self):
      process = self.get_process()
      #logger.info('process = '+ process)
      if process is None:
         return None
      njets = None
      for par in self.parameters:
         if 'njets' in par.comment:
            njets = int(par.value)
            break
      #logger.info('njets = '+ str(njets))
      if njets is not None:
         return AlpgenParameters.event_sizes[process][njets]
      return None

   def __str__(self):
      if self.filename is None:
         s = ' no file set '
         return s
      s = ''
      s += 'Alpgen Unweighted Parameters File: \n'
      s += '   ' + self.filename + '\n'
      s += ' hard process code = %3i (process = %s)\n' % (self.hard_process_code,self.get_process())
      s += ' m_c = %5.2f m_b = %5.2f m_t = %5.2f\n' % (self.mc,self.mb,self.mt)
      s += ' m_w = %5.2f m_z = %5.2f m_h = %5.2f\n' % (self.mw,self.mz,self.mh)
      s += ' cross section =  %10.8e +/- %10.8e\n' % (self.cross_section,self.cross_section_error)
      s += ' luminosity = %10.8e  event count = %15i\n' % (self.lumi,self.event_count)
      for par in self.parameters:
         s += str(par) + '\n'
      return s

   @staticmethod
   def read_file(filename):
      if not os.path.exists(filename):
         logger.error('cannot read file ' + filename)
         return None

      par = AlpgenUnwParFile()
      par.filename = filename
      in_par = False
      before_par = True
      for line in open(filename):
         if AlpgenUnwParFile.BEGIN_PARAMETERS in line:
            in_par = True
            before_par = False
            continue
         elif AlpgenUnwParFile.END_PARAMETERS in line:
            in_par = False
            continue
         elif AlpgenUnwParFile.HARD_PROCESS_ID in line:
            par.hard_process_code = int(line.split()[0])
            continue
         elif AlpgenUnwParFile.MASSES_ID in line:
            split_str = line.split()
            par.mc = float(split_str[0])
            par.mb = float(split_str[1])
            par.mt = float(split_str[2])
            par.mw = float(split_str[3])
            par.mz = float(split_str[4])
            par.mh = float(split_str[5])
            continue
         elif AlpgenUnwParFile.CROSS_SECTION_ID in line:
            split_str = line.split()
            par.cross_section = float(split_str[0])
            par.cross_section_error = float(split_str[1])
            continue
         elif AlpgenUnwParFile.EVENT_COUNT_ID in line:
            split_str = line.split()
            par.event_count = int(split_str[0])
            par.lumi = float(split_str[1])
            continue

         if before_par:
            par.start_lines.append(line)

         if in_par:
            par.parameters.append(AlpgenUnwParFile.AlpgenUnwParameter.read_line(line[0:-1]))

      return par

   def add_pythia6_settings(self,jet_treatment=EXCLUSIVE_JET_MATCHING):#,output_filename=None):
      
      #if output_filename is None:
      #   output_filename = self.filename + '.tmp'

      new_parameters_list = [
                  AlpgenUnwParFile.AlpgenUnwParameter(501,20.,'min ETCLUS used for parton-jet matching (Normally ETCLUS = ptjmin + 5)'),
                  AlpgenUnwParFile.AlpgenUnwParameter(502,self.get_drjmin(),'min RCLUS value for parton-jet matching'),
                  AlpgenUnwParFile.AlpgenUnwParameter(503,6.0,'max ETACLUS value for parton-jet matching'),
                  AlpgenUnwParFile.AlpgenUnwParameter(504,jet_treatment,'Jet Matching: 0 inclusive 1 exclusive'),
                 ]

      # check that these parameters do not already exist
      for parameter in self.parameters:
         for new_par in new_parameters_list:
            if new_par.key == parameter.key:
               logger.error('Key ' + str(new_par.key) + ' already exits in current parameter file. Stopping update.')
               raise Exception('Key ' + str(new_par.key) + ' already exits in current parameter file. Stopping update.')

      self.parameters += new_parameters_list
      #self.write_file(output_filename)
      #try:
      #   shutil.copy(self.filename,input_filename+'.old')
      #   os.remove(self.filename)
      #   shutil.copy(output_filename,self.filename)
      #   os.remove(output_filename)
      #except:
      #   logger.error(' exception caught: ' + str(sys.exc_info()[1]))
      #   raise


if __name__ == '__main__':
   import optparse
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='')
   parser.add_option('-i','--input-file',dest='input_filename',help='Input Filename to test')
   options,args = parser.parse_args()

   
   manditory_args = [
                     'input_filename',
                  ]

   for man in manditory_args:
      if not options.__dict__[man]:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)

   if not os.path.exists(options.input_filename):
      logger.error('file does not exist: ' + options.input_filename)
      sys.exit(-1)

   parfile = AlpgenUnwParFile.read_file(options.input_filename)
   logger.info('\n' + str(parfile) + '\n')



