#!/usr/bin/env python
import os,sys,optparse,logging,glob,filecmp
from sherpa.SherpaInputCard import SherpaInputCard
logger = logging.getLogger(__name__)

def main():
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='Make sure a sherpa card has the needed settings.')
   parser.add_option('-i','--input',dest='input',help='Input sherpa card to modify.')
   parser.add_option('-g','--glob',dest='input_pattern',help='A glob that will be used to get a list of input cards for modification. If this option is used, "-i" is ignored and the "-o" option is used as a new folder where the modified files will be written.')
   parser.add_option('-o','--output',dest='output',help='output sherpa card')
   parser.add_option('-n','--nevts',dest='nevts',help='Number of events (per rank if running MPI) to include in the card.',type='int')
   options,args = parser.parse_args()

   
   manditory_args = [
                     'output',
                     'nevts',
                  ]

   for man in manditory_args:
      if not options.__dict__[man]:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)

   if ( (options.input is None and options.input_pattern is None) or
        (options.input is not None and options.input_pattern is not None)):
      logger.error('Must specify either input or glob pattern, not both or none.')
      parser.print_help()
      sys.exit(-1)
   
   
   if options.input is not None:
      modify_card(options.input,options.output,options.nevts)
   else:
      input_cards = glob.glob(options.input_pattern)

      for input_card in input_cards:
         output = os.path.join(options.output,os.path.basename(input_card))
         # make sure card doesn't exist
         if not os.path.exists(output):
            modify_card(input_card,output,options.nevts)
         else:
            logger.warning(' file ' + output + ' already exists, not overwriting.')


def modify_card(input,output,nevts):
   input_card = SherpaInputCard.read_file(input)

   # add lines to input card in (run) block
   run_block = input_card.blocks[SherpaInputCard.RUN_BLOCK]
   run_block.add_line('')
   if not run_block.contains('EVENTS'):
      run_block.add_line('EVENTS %i;' % nevts)
   if not run_block.contains('EVENT_OUTPUT'):
      run_block.add_line('EVENT_OUTPUT=HepMC_GenEvent[integrate];')
   run_block.add_line('')
   if not run_block.contains('% collider setup'):
      run_block.add_line('% collider setup')
   if not run_block.contains('BEAM_1'):
      run_block.add_line('BEAM_1 2212; BEAM_2 2212;')
   if not run_block.contains('BEAM_ENERGY_1'):
      run_block.add_line('BEAM_ENERGY_1 = 6500.; BEAM_ENERGY_2 = 6500.;')
   run_block.add_line('')
   if not run_block.contains('% model parameters'):
      run_block.add_line('% model parameters')
   if not run_block.contains('MASS[6]'):
      run_block.add_line('MASS[6]=172.5;')
   if not run_block.contains('MASS[23]'):
      run_block.add_line('MASS[23]=91.1876;')
   if not run_block.contains('MASS[24]'):
      run_block.add_line('MASS[24]=80.399;')
   if not run_block.contains('WIDTH[23]'):
      run_block.add_line('WIDTH[23]=2.4952;')
   if not run_block.contains('WIDTH[24]'):
      run_block.add_line('WIDTH[24]=2.085;')
   if not run_block.contains('SIN2THETAW'):
      run_block.add_line('SIN2THETAW=0.23113;')
   if not run_block.contains('SIN2THETAW'):
      run_block.add_line('SIN2THETAW=0.23113;')
   if not run_block.contains('MAX_PROPER_LIFETIME'):
      run_block.add_line('MAX_PROPER_LIFETIME=10.0;')
   if not run_block.contains('WIDTH[15]'):
      run_block.add_line('WIDTH[15]=2.26735e-12;')

   # remove block (model)
   input_card.delete_block(SherpaInputCard.MODEL_BLOCK)

   input_card.write_file(output)

if __name__ == "__main__":
   main()
