#!/usr/bin/env python
import os,sys,subprocess,logging,optparse,shutil
logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')
logger = logging.getLogger(__name__)
# import tools
sys.path.append('/users/hpcusers/svn/tools/python/athena')
from alpgen_unw_dataset import create_alpgen_dataset
from prun_gen_trf_grid import prun

TMP_PRUN_FOLDER='tmp_prun'
GRID_SCRIPT='/users/hpcusers/svn/tools/python/athena/gen_trf_grid.py'

def main():
   # print command
   print str(sys.argv)
   # parse command line
   parser = optparse.OptionParser(description='finish a job on the grid with Generate_trf.py')
   parser.add_option('-i','--input-base',dest='input_basename',help='input base filename for the alpgen ".unw" and "_unw.par" files, can be a comma separated list of base names with pathes (alpout)',default='alpout')
   parser.add_option('-d','--inDS',dest='inDS',help='name of the dataset which will be created from the alpgen data')
   parser.add_option('-q','--outDS',dest='outDS',help='Output dataset name ("")',default='')
   parser.add_option('-k','--nEvtsPerFile',dest='nEvtsPerFile',help='if splitting the input alpgen file, how many events per file.(5000)',type='int',default=5000)
   parser.add_option('-a','--asetup',dest='asetup',help='athena asetup configuration string (17.2.11.12,AltasProduction)',default='17.2.11.12,AtlasProduction')
   parser.add_option('-s','--site',dest='site',help='Host site for the uploaded dataset (MWT2_UC_SCRATCHDISK)',default='MWT2_UC_SCRATCHDISK')
   parser.add_option('-e','--ecmEnergy',dest='ecmEnergy',help='Center of Mass Energy for this sample (8000)',default='8000')
   parser.add_option('-r','--randomSeed',dest='randomSeed',help='Random Seed to use (12345678)',type='int',default=12345678)
   parser.add_option('-o','--outputEVNTFile',dest='outputEVNTFile',help='output EVNT filename (myOutputEVNTFile.root)',default='myOutputEVNTFile.root')
   parser.add_option('-f','--firstEvent',dest='firstEvent',help='First event to start on (1)',type='int',default=1)
   parser.add_option('-n','--runNumber',dest='runNumber',help='Run Number')
   parser.add_option('-c','--jobConfig',dest='jobConfig',help='Job configuration, AKA Starting Job Options')
   parser.add_option('-j','--evgenJobOpts',dest='evgenJobOpts',help='Input tarball with all the generator control job options')
   parser.add_option('-x','--skip-dataset',dest='createDataset',help='Use this flag to skip the dataset creation.',action='store_false',default='True')
   
   options,args = parser.parse_args()

   if options.inDS is None:
      parser.error('-d is required')
   if options.runNumber is None:
      parser.error('-n is required')
   if options.jobConfig is None:
      parser.error('-c is required')
   if options.evgenJobOpts is None:
      parser.error('-j is required')

   if options.outDS == '':
      options.outDS = options.inDS + '.gentrf'

   logger.info(' Running with settings: ')
   logger.info('    input base         = ' + str(options.input_basename.split(',')))
   logger.info('    n evts per file    = ' + str(options.nEvtsPerFile))
   logger.info('    site               = ' + options.site)
   logger.info('    asetup             = ' + options.asetup)
   logger.info('    ecmEnergy          = ' + options.ecmEnergy)
   logger.info('    randomSeed         = ' + str(options.randomSeed))
   logger.info('    outputEVNTFile     = ' + options.outputEVNTFile)
   logger.info('    runNumber          = ' + options.runNumber)
   logger.info('    firstEvent         = ' + str(options.firstEvent))
   logger.info('    jobConfig          = ' + options.jobConfig)
   logger.info('    evgenJobOpts       = ' + options.evgenJobOpts)
   logger.info('    inDS               = ' + options.inDS)
   logger.info('    outDS              = ' + options.outDS)

   # this function splits the alpgen input files and uploads the tarball as a data set
   if options.createDataset:
      logger.info(' running create_alpgen_dataset ')
      try:
         create_alpgen_dataset(options.input_basename.split(','),options.nEvtsPerFile,options.inDS,options.site,options.asetup)
      except:
         logger.error(' exception caught: ' + str(sys.exc_info()[1]))
         return
   else:
      logger.info(' skipping create_alpgen_dataset ')

   logger.info(' running prun ')
   # run a prun job which runs Generate_trf.py
   try:
      prun(
            options.asetup,
            options.ecmEnergy,
            options.randomSeed,
            options.outputEVNTFile,
            options.runNumber,
            options.firstEvent,
            options.jobConfig,
            options.evgenJobOpts,
            options.inDS,
            options.outDS,
            GRID_SCRIPT,
          )
   except:
      logger.error(' exception caught: ' + str(sys.exc_info()[1]))
      return

if __name__ == '__main__':
   main()
