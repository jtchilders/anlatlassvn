#!/usr/bin/python
import os,sys,optparse,logging
logging.basicConfig(filename='out.log',format='%(asctime)s %(levelname)s:%(name)s:%(message)s',level=logging.INFO)
logger = logging.getLogger(__name__)

XRDCP='xrdcp'
GENERATE_TRF='Generate_trf.py'

def main():
   parser = optparse.OptionParser(description='run Generate_trf.py for Alpgen+Pythia job')
   parser.add_option('-e','--ecmEnergy',dest='ecmEnergy',help='Center of Mass Energy for this sample')
   parser.add_option('-r','--randomSeed',dest='randomSeed',help='Random Seed to use',type='int',default=12345678)
   parser.add_option('-o','--outputEVNTFile',dest='outputEVNTFile',help='output EVNT filename',default='output.root')
   parser.add_option('-f','--firstEvent',dest='firstEvent',help='First event to start on',type='int',default=1)
   parser.add_option('-n','--runNumber',dest='runNumber',help='Run Number')
   parser.add_option('-c','--jobConfig',dest='jobConfig',help='Job configuration, AKA Starting Job Options')
   parser.add_option('-i','--inputGeneratorFile',dest='inputGeneratorFile',help='Input tarball with generated events and parameters')
   parser.add_option('-j','--evgenJobOpts',dest='evgenJobOpts',help='Input tarball with all the generator control job options')
   options,args = parser.parse_args()

   if options.ecmEnergy is None:
      parser.error('Must set -e')
   if options.runNumber is None:
      parser.error('Must set -n')
   if options.jobConfig is None:
      parser.error('Must set -c')
   if options.inputGeneratorFile is None:
      parser.error('Must set -i')
   if options.evgenJobOpts is None:
      parser.error('Must set -j')


   logger.info(' Running with settings: ')
   logger.info('    ecmEnergy          = ' + options.ecmEnergy)
   logger.info('    randomSeed         = ' + str(options.randomSeed))
   logger.info('    outputEVNTFile     = ' + options.outputEVNTFile)
   logger.info('    runNumber          = ' + options.runNumber)
   logger.info('    firstEvent         = ' + str(options.firstEvent))
   logger.info('    jobConfig          = ' + options.jobConfig)
   logger.info('    inputGeneratorFile = ' + options.inputGeneratorFile)
   logger.info('    evgenJobOpts       = ' + options.evgenJobOpts)
   
   # copy dataset file to local work area
   cmd = XRDCP + ' ' + options.inputGeneratorFile + ' ./' + os.path.basename(options.inputGeneratorFile)
   logger.info(' copy data file to local area ' + cmd)
   os.system(cmd)
   # untar evgenJobOpts
   cmd = 'tar zxf ' + options.evgenJobOpts
   logger.info(' untar ' + options.evgenJobOpts)
   os.system(cmd)
   cmd = ( GENERATE_TRF + 
            ' ecmEnergy=' + options.ecmEnergy + 
            ' randomSeed=' + str(options.randomSeed) + 
            ' outputEVNTFile=' + options.outputEVNTFile + 
            ' runNumber=' + options.runNumber + 
            ' firstEvent=' + str(options.firstEvent) + 
            ' jobConfig=' + options.jobConfig + 
            ' inputGeneratorFile=' + os.path.basename(options.inputGeneratorFile)
          )
   logger.info('command: ' + cmd)
   os.system(cmd)
   
   #os.system('tar zcf ' + options.outputEVNTFile + '.tgz ' + options.outputEVNTFile)

   return 0


if __name__ == "__main__":
   sys.exit(main())




