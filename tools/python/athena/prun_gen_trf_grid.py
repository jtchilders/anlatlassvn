#!/usr/bin/python
import os,sys,optparse,logging,subprocess,shlex,shutil
logging.basicConfig(format='%(asctime)s %(levelname)s:%(name)s:%(message)s',level=logging.INFO)
logger = logging.getLogger(__name__)

SCRIPT_ON_GRID='gen_trf_grid.py'
TMP_FOLDER='tmp_prun'

def main():
   parser = optparse.OptionParser(description='submit grid job with prun for Alpgen+Pythia job')
   parser.add_option('-a','--asetup',dest='asetup',help='config string for athena asetup command')
   parser.add_option('-e','--ecmEnergy',dest='ecmEnergy',help='Center of Mass Energy for this sample')
   parser.add_option('-r','--randomSeed',dest='randomSeed',help='Random Seed to use',type='int',default=12345678)
   parser.add_option('-o','--outputEVNTFile',dest='outputEVNTFile',help='output EVNT filename',default='myOutputEVNTFile.root')
   parser.add_option('-f','--firstEvent',dest='firstEvent',help='First event to start on',type='int',default=1)
   parser.add_option('-n','--runNumber',dest='runNumber',help='Run Number')
   parser.add_option('-c','--jobConfig',dest='jobConfig',help='Job configuration, AKA Starting Job Options')
   parser.add_option('-j','--evgenJobOpts',dest='evgenJobOpts',help='Input tarball with all the generator control job options')
   parser.add_option('-p','--inDS',dest='inDS',help='Input dataset name')
   parser.add_option('-q','--outDS',dest='outDS',help='Output dataset name')
   options,args = parser.parse_args()
   
   if options.asetup is None:
      parser.error('Must set -a')
   if options.ecmEnergy is None:
      parser.error('Must set -e')
   if options.runNumber is None:
      parser.error('Must set -n')
   if options.jobConfig is None:
      parser.error('Must set -c')
   if options.evgenJobOpts is None:
      parser.error('Must set -j')
   if options.inDS is None:
      parser.error('Must set -p')
   if options.outDS is None:
      parser.error('Must set -q')

   logger.info(' Running with settings: ')
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
   
   prun(options.asetup,
        options.ecmEnergy,
        options.randomSeed,
        options.outputEVNTFile,
        options.runNumber,
        options.firstEvent,
        options.jobConfig,
        options.evgenJobOpts,
        options.inDS,
        options.outDS
       )

def prun(
         asetup,
         ecmEnergy,
         randomSeed,
         outputEVNTFile,
         runNumber,
         firstEvent,
         jobConfig,
         evgenJobOpts,
         inDS,
         outDS,
         grid_script=SCRIPT_ON_GRID,
         ATLAS_LOCAL_ROOT_BASE = '/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase',
         TMP_SHELL_SCRIPT = 'tmp_prun_gen_trf_grid.sh'
        ):
   
   # get current working directory
   starting_folder = os.getcwd()

   # create a local temporary working area to seutp the job
   # create a folder in which to run prun
   if os.path.exists(TMP_FOLDER):
      logger.info(' Temporary PRUN folder, ' + TMP_FOLDER + ', already exists')
   else:
      os.mkdir(TMP_FOLDER)
   # copy the needed files to that folder
   if os.path.exists(jobConfig):
      shutil.copy(jobConfig,TMP_FOLDER+'/')
      # strip path as it is no longer needed
      jobConfig = os.path.basename(jobConfig)
   else:
      logger.info(' jobConfig: ' + jobConfig + ' does not exist locally, will try to continue ')
   if os.path.exists(evgenJobOpts):
      shutil.copy(evgenJobOpts,TMP_FOLDER+'/')
      # strip path as it is no longer needed
      evgenJobOpts = os.path.basename(evgenJobOpts)
   else:
      logger.info(' evgenJobOpts: ' + evgenJobOpts + ' does not exist locally, will try to continue ')
   shutil.copy(grid_script,TMP_FOLDER+'/')
   # move to working directory
   os.chdir(TMP_FOLDER)

   if not os.path.exists(os.path.basename(grid_script)):
      logger.error('The script that is run on the grid, ' + os.path.basenameg(grid_script) + ', is not in current directory')
      raise Exception(' grid script is missing ')

   # first need to setup Athena
   out = ("""#!/usr/bin/env bash
# setup ASetup environment
export ATLAS_LOCAL_ROOT_BASE=""" + ATLAS_LOCAL_ROOT_BASE + """
# setup ASetup environment
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
echo ' setting up prun '
localSetupPandaClient
# run prun command
echo ' running prun command '
prun --exec="python """ + os.path.basename(grid_script) + ' -e ' + str(ecmEnergy) + ' -o ' + outputEVNTFile + """ \\
     -n """ + str(runNumber) + ' -c ' + jobConfig + """ \\
     -i %IN -j """ + evgenJobOpts + """ " \\
     --inDS=""" + inDS + """ \\
     --nFilesPerJob=1 \\
     --outDS=""" + outDS + """ \\
     --outputs=""" + outputEVNTFile + """ \\
     --athenaTag=""" + asetup + '\n'
   )
   
   logger.debug(' file contents: \n' + out )

   file = open(TMP_SHELL_SCRIPT,'w')
   file.write(out)
   file.close()
   file = None
   
   logger.info(' launching prun script ')
   p = subprocess.Popen('sh ' + TMP_SHELL_SCRIPT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
   line = p.stdout.readline()
   while line != '':
      logger.info(line[0:-1])
      line = p.stdout.readline()
   p.wait()
   if p.returncode != 0:
      stdout,stderr = p.communicate()
      logger.error('ERROR Setting up athena and running prun command. returncode = ' + str(p.returncode) + ' stderr:\n' + stderr)
      # return to working directory
      os.chdir(starting_folder)
      raise Exception(' prun command failed ')
   # return to working directory
   os.chdir(starting_folder)
   logger.debug(' finished prun ' )



if __name__ == "__main__":
   sys.exit(main())


