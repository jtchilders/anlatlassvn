#!/usr/bin/env python
import sys,os,logging,optparse
logger = logging.getLogger(__name__)


class SherpaInputCard:
   RUN_BLOCK         = '(run)'
   PROCESSES_BLOCK   = '(processes)'
   SELECTOR_BLOCK    = '(selector)'
   MODEL_BLOCK       = '(model)'
   BLOCK_BEGIN       = '{'
   BLOCK_END         = '}'
   BLOCKS            = [ RUN_BLOCK,
                         PROCESSES_BLOCK,
                         SELECTOR_BLOCK,
                         MODEL_BLOCK,
                       ]

   RANDOM_SEED       = 'RANDOM_SEED'
   EVENTS            = 'EVENTS'

   class Setting:
      def __init__(self,parameter = None, value = None, operator = None):
         self.parameter = parameter
         self.value     = value
         self.operator  = operator
      def fill(self,string):
         string = string.strip()
         #logger.info('filling string: ' + string)
         if ':=' in string:
            split = string.split(':=')
            self.parameter    = split[0].strip()
            self.value        = split[1].strip()
            self.operator     = ':='
         elif '=' in string:
            split = string.split('=')
            self.parameter    = split[0].strip()
            self.value        = split[1].strip()
            self.operator     = '='
         else:
            split = string.split(' ')
            self.parameter    = split[0].strip()
            self.value        = ''
            for i in range(len(split)-1):
               self.value     += split[i+1]
               if i < len(split)-2:
                  self.value     += ' '
            self.operator     = None
         #logger.info('    got: ' + self.parameter + ' ' + str(self.operator) + ' ' + self.value)
      def __str__(self):
         if self.operator is not None:
            return self.parameter + self.operator + self.value
         return self.parameter + ' ' + self.value
      def get_parameter_index(self,parameter_name):
         index = -1
         for setting in self.settings:
            if type(setting) is str:
               pass
            elif isinstance(setting,SherpaInputCard.Setting):
               if parameter_name == setting.parameter:
                  return self.settings.index(setting)
         return index


   class Block:
      def __init__(self,blockname):
         self.settings = []
         self.blockname = blockname
      def add_line(self,line):
         line = line.replace('\n','')
         line = line.strip()
         if len(line) > 0:
            if ';' in line:
               split = line.split(';')
               for setting_string in split:
                  if len(setting_string.lstrip()) > 0:
                     setting = SherpaInputCard.Setting()
                     setting.fill(setting_string)
                     self.settings.append(setting)
            elif line[0] == '%':
               self.settings.append(line)
            else:
               self.settings.append(line)
         
      def __str__(self):
         s = ''
         s += '%s%s\n' % (self.blockname,SherpaInputCard.BLOCK_BEGIN)
         for setting in self.settings:
            if type(setting) is str: # setting is a comment
               s += '\n  %s\n' % setting
            elif isinstance(setting,SherpaInputCard.Setting):
               s += '  %s;\n' % str(setting)
         s += '%s%s\n' % (SherpaInputCard.BLOCK_END,self.blockname)
         return s


   def __init__(self):
      self.blocks = {}
      pass

   @staticmethod 
   def read_file(filename):
      card = None
      if os.path.exists(filename):
         # create empty card
         card = SherpaInputCard()
         # need to know which block we are in
         current_block = None
         
         # open input file
         for line in open(filename):
            in_transition = False
            # check for block transition
            for block in SherpaInputCard.BLOCKS:
               if block in line:
                  in_transition = True
                  if SherpaInputCard.BLOCK_BEGIN in line:
                     current_block = block
                     card.blocks[block] = SherpaInputCard.Block(block)
                  elif SherpaInputCard.BLOCK_END in line:
                     current_block = None
            # if the current block is sent, save this line to that block
            if current_block is not None and not in_transition:
               card.blocks[current_block].add_line(line.replace('\n','').lstrip())
      else:
         logger.error('No File Found: ' + str(filename))

      return card

   def write_file(self,filename):
      outfile = open(filename,'w')
      for block in self.blocks.values():
         outfile.write(str(block))
         outfile.write('\n')

   def __str__(self):
      s = ''
      for block in self.blocks.values():
         s += str(block)
         s += '\n'
      return s

   def delete_block(self,blockname):
      if blockname in self.blocks.keys():
         del self.blocks[blockname]

   def set_random_seeds(self,seed1,seed2,seed3,seed4):
      run_block = self.blocks[SherpaInputCard.RUN_BLOCK]
      index = run_block.get_parameter_index(SherpaInputCard.RANDOM_SEED)
      seed_string = str(seed1) + ' ' + str(seed2) + ' ' + str(seed3) + ' ' + str(seed4)
      if index < 0: # no random seed set
         run_block.settings.append(SherpaInputCard.Setting(SherpaInputCard.RANDOM_SEED,'=',seed_string))
      else:
         run_block.settings[index].value = seed_string

   def set_number_events(self,number_events):
      run_block = self.blocks[SherpaInputCard.RUN_BLOCK]
      index = run_block.get_parameter_index(SherpaInputCard.EVENTS)
      if index < 0:
         run_block.settings.append(SherpaInputCard.Setting(SherpaInputCard.EVENTS,'=',str(number_events)))
      else:
         run_block.settings[index].value = str(number_events)



if __name__ == '__main__':
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='test card copying')
   parser.add_option('-i','--input',dest='input',help='input card')
   parser.add_option('-o','--output',dest='output',help='output copy of input')
   options,args = parser.parse_args()

   
   manditory_args = [
                     'input',
                     'output',
                  ]

   for man in manditory_args:
      if options.__dict__[man] is None:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)

   card = SherpaInputCard.read_file(options.input)
   print str(card)
   card.write_file(options.output)

