import logging
import re
logger = logging.getLogger(__name__)

class InputFile:
   ''' Class for configuring application and can write a cmd file '''
      
   def __init__(self):
      # Mandatory parameters:
      self.filename_base   = 'sherpa'
      self.nevt            = 100000 # number of events per iteration

      # Options for generation
      self.options = []

   def __str__(self):
      txtfmt = '''
(run){
    EVENTS %d; ERROR 0.99;
    EVENT_OUTPUT=HepMC_Short[%s];

    %% scales, tags for scale variations
    SP_NLOCT 1; FSF:=1.; RSF:=1.; QSF:=1.;
    SCALES METS{FSF*MU_F2}{RSF*MU_R2}{QSF*MU_Q2};

    %% tags for process setup
    NJET:=2; LJET:=2; QCUT:=20.;

    %% me generator settings
    ME_SIGNAL_GENERATOR Comix Amegic LOOPGEN;
    EVENT_GENERATION_MODE Weighted;
    LOOPGEN:=Internal;

    %% exclude tau from lepton container
    MASSIVE[15] 1;

    %% collider setup
    BEAM_1 2212; BEAM_ENERGY_1 = 4000.;
    BEAM_2 2212; BEAM_ENERGY_2 = 4000.;

    AMEGIC_LIBRARY_MODE 1; 
}(run)

(processes){
    Process 93 93 -> 90 91 93{NJET};
    Order_EW 2; CKKW sqr(QCUT/E_CMS);
    NLO_QCD_Mode MC@NLO {LJET};
    ME_Generator Amegic {LJET};
    Loop_Generator LOOPGEN {LJET};
    Scales LOOSE_METS{FSF*MU_F2}{RSF*MU_R2}{QSF*MU_Q2} {7,8};
    End process;
}(processes)

(selector){
    Mass 11 -12 5. E_CMS
    Mass 13 -14 5. E_CMS
    Mass -11 12 5. E_CMS
    Mass -13 14 5. E_CMS
}(selector)
'''
      txt = txtfmt % (self.nevt,self.filename_base)

      return txt

   def write(self,filename):
      try:
         fp = open(filename,'w')
      except IOError,e:
         logger.error('opening filename '+filename+' for writing. Exception: '+str(e))
         return -1
      
      try:
         fp.write(str(self))
      except IOError,e:
         logger.error('writing to filename ' + filename + ' for writing. Exception: ' + str(e))
         return -2

      fp.close()
      return 0

   def edit_option(self,label,value):
      for option in self.options:
         if option.label == label:
            option.value = value
            break

   def read(self,filename):
      try:
         fp = open(filename)
      except IOError,e:
         logger.error('opening filename '+filename+' for reading. Exception: ' + str(e))
         return -3

      filelines = fp.readlines()
      filelines = [f.strip() for f in filelines]
      events_lines = [f for f in filelines if f.startswith('EVENTS')]
      if events_lines:
         self.nevt = events_lines[0].split()[1]
      event_output_lines = [f for f in filelines if f.startswith('EVENT_OUTPUT')]
      if event_output_lines:
         result = re.match('.*[(\d+)]',event_output_lines[0])
         print 'result = ', result
         if result:
            self.filename_base = result.group(1)
      return 0

if __name__ == '__main__':
   import sys
   f = InputFile()
   f.nevt = 999
   f.filename_base = 'sherpa_test'
   if len(sys.argv)>1:
     f.read(sys.argv[1])
   print f.filename_base, f.nevt, f


         
