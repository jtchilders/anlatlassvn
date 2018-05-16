#!/usr/bin/env python
import os,sys,optparse,logging,subprocess
logger = logging.getLogger(__name__)
from athena.run_athena_in_shell import run_athena_in_shell

ASETUP_ARGS="--testarea=/users/hpcusers/athenaInstallArea 17.7.3.6,AtlasProduction,gcc46"
EXECUTABLE='Generate_trf.py'
ARGS='ecmEnergy=%i randomSeed=%i outputEVNTFile=%s runNumber=%i firstEvent=%i jobConfig=%s inputGeneratorFile=%s evgenJobOpts=%s postExec="topAlg.CountHepMC.RequestedOutput=%i"'
DEFAULT_RANDOM_SEED = 131376
DEFAULT_FIRST_EVENT = 1
DEFAULT_REQUESTED_OUTPUT = 100000

def main():
   logging.basicConfig(level=logging.INFO)
   parser = optparse.OptionParser(description='Run Generate_trf.py over an alpgen data input file.')
   parser.add_option('-e','--ecmEnergy',dest='ecmEnergy',help='Center-of-mass energy.',type=int)
   parser.add_option('-o','--outputEVNTFile',dest='outputEVNTFile',help='Output EVNT file name.')
   parser.add_option('-r','--runNumber',dest='runNumber',help='Run Number.',type=int)
   parser.add_option('-s','--randomSeed',dest='randomSeed',help='Random seed [DEFAULT = ' + str(DEFAULT_RANDOM_SEED) + '] ',type=int,default=DEFAULT_RANDOM_SEED)
   parser.add_option('-f','--firstEvent',dest='firstEvent',help='First event at which to start [DEFAULT = ' + str(DEFAULT_FIRST_EVENT) + ']',type=int,default=DEFAULT_FIRST_EVENT)
   parser.add_option('-c','--jobConfig',dest='jobConfig',help='Job Config (job options file).')
   parser.add_option('-i','--inputGeneratorFile',dest='inputGeneratorFile',help='Input Generator File from Alpgen. This should be a tarball that follows the naming pattern "*alpgen.XXXXXX.TXT.v1*.tgz" which contains two files, one should be named as "*alpgen.XXXXXX.TXT.v1*.dat" (this is the _unw.par file) and one should be named as "*alpgen.XXXXXX.TXT.v1*.events" (this is the .unw file).')
   parser.add_option('-j','--evgenJobOpts',dest='evgenJobOpts',help='Event Generation Job Options Tarball name.')
   parser.add_option('-n','--minevents',dest='minevents',help='Generate_trf requires you to tell it the minimimum events you want. In this case we just set this value too high such that all input events are processed. The unfortunate side effect of this is that Athena exits with an error since it does not get the minimum number of events. [DEFAULT = ' + str(DEFAULT_REQUESTED_OUTPUT) + '] ',type=int,default=DEFAULT_REQUESTED_OUTPUT)
   
   options,args = parser.parse_args()
 
   logger.info('==============================================================')
   logger.info(' running with following options: ')
   logger.info(' ecmEnergy          = ' + str(options.ecmEnergy) )
   logger.info(' outputEVNTFile     = ' + str(options.outputEVNTFile) )
   logger.info(' runNumber          = ' + str(options.runNumber) )
   logger.info(' randomSeed         = ' + str(options.randomSeed) )
   logger.info(' firstEvent         = ' + str(options.firstEvent) )
   logger.info(' jobConfig          = ' + str(options.jobConfig) )
   logger.info(' inputGeneratorFile = ' + str(options.inputGeneratorFile) )
   logger.info(' evgenJobOpts       = ' + str(options.evgenJobOpts) )
   logger.info(' minevents          = ' + str(options.minevents) )
   logger.info('==============================================================')

   if options.ecmEnergy is None:
      parser.error('must specify center of mass energy.')
   if options.outputEVNTFile is None:
      parser.error('must specify output EVNT file name.')
   if options.jobConfig is None:
      parser.error('must specify job config file.')
   if options.runNumber is None:
      parser.error('must specify run number.')
   if options.evgenJobOpts is None:
      options.evgenJobOpts = ''
      #parser.error('must specify event generation job options.')
   if options.inputGeneratorFile is None:
      parser.error('must specify input generator file.')
  

   return run_athena_in_shell( 
                        EXECUTABLE,
                        ARGS % (options.ecmEnergy,options.randomSeed,options.outputEVNTFile,
                                options.runNumber,options.firstEvent,options.jobConfig,
                                options.inputGeneratorFile,options.evgenJobOpts,options.minevents
                               ),
                        ASETUP_ARGS,
                       )

if __name__ == '__main__':
   sys.exit(main())


