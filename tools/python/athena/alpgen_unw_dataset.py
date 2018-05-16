#!/usr/bin/python
import os,sys,optparse,logging,subprocess,shlex
logging.basicConfig(format='%(asctime)s %(levelname)s:%(name)s:%(message)s',level=logging.INFO)
logger = logging.getLogger(__name__)

TMP_FOLDER='tmp_alpgen_unw_split'
TMP_SPLIT_FILEBASE='split.'
TARBALL_FOLDER='dataset'
TMP_SHELL_SCRIPT='tmp_upload_dataset.sh'

from AlpgenUnwFile import AlpgenUnwFile,ALPGEN_ATHENA_FILEBASE
from AlpgenUnwParFile import AlpgenUnwParFile

def main():
   parser = optparse.OptionParser(description='divide alpgen unweighted events file and parameters file, create tarball for each subset and upload to the grid as a dataset')
   parser.add_option('-a','--asetup',dest='asetup',help='athena asetup config string')
   parser.add_option('-i','--input-base',dest='inbase',help='base filename for the input alpgen files, can be a comma separated list',default='alpout')
   parser.add_option('-n','--nevts',dest='nevts',help='number of events per sub-file',type='int',default=5000)
   parser.add_option('-o','--outDS',dest='outDS',help='output dataset name which is uploaded to the grid')

   parser.add_option('-s','--site',dest='site',help='site to upload the dataset',default='MWT2_UC_SCRATCHDISK')

   
   options,args = parser.parse_args()
   
   if options.outDS is None:
      parser.error('must set -o')
   if options.asetup is None:
      parser.error('must set -a')

   logger.info(' Running with settings: ')
   logger.info('    asetup             = ' + options.asetup)
   logger.info('    input base         = ' + options.inbase)
   logger.info('    events per file    = ' + str(options.nevts))
   logger.info('    output DS          = ' + options.outDS)
   logger.info('    site               = ' + options.site)
   
   create_alpgen_dataset(
        options.inbase.split(','),
        options.nevts,
        options.outDS,
        options.site,
        options.asetup,
       )


''' 
   This function takes Alpgen unweighted event files and splits them into smaller
   files of 'nevts' each. These files are then uploaded to the grid as dataset
   'outDS' and stored at 'site'.
'''

def create_alpgen_dataset(
         input_basenames, #  list []
         nevts, # int
         outDS, # string
         site, # string
         asetup, # string
         ATLAS_LOCAL_ROOT_BASE='/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase',
        ):

   # make a temporary directory to hold all the files
   try:
      os.mkdir(TMP_FOLDER)
   except:
      logger.error('Error trying to create folder ' + TMP_FOLDER + ' Exception: ' + str(sys.exc_info()[1]) + '. Will Try to continue, but it may fail')

   # first divide the alpgen unweighted events/parameter file into smaller bits
   split_basenames = []
   for input_basename in input_basenames:
      logger.info(' splitting files with basename: ' + input_basename)
      alpgenUnwParFile = AlpgenUnwParFile(input_basename+'_unw.par')
      alpgenUnwFile = AlpgenUnwFile(input_basename+'.unw',alpgenUnwParFile)
      # make a tmp folder for the files
      try:
         split_basenames += alpgenUnwFile.SplitFile(nevts,
                                                    os.path.join(TMP_FOLDER,TMP_SPLIT_FILEBASE),
                                                    file_number_offset=len(split_basenames),
                                                   )
      except:
         logger.error('Error while splitting alpgen file ' + alpgenUnwFile.filename + ' Exception: ' + str(sys.exc_info()[1]) )
         raise
         
   # now make the tarball of each unw/unw.par file set
   logger.info(' looping over all files and making tarball of each .unw and _unw.par set. ')
   try:
      os.mkdir(os.path.join(TMP_FOLDER,TARBALL_FOLDER))
      alpgenUnwFile.MakeDataSetTarball(
            split_basenames,
            os.path.join(TMP_FOLDER,TARBALL_FOLDER,ALPGEN_ATHENA_FILEBASE+outDS+'.')
          )
   except:
      logger.error('Error trying to create folder ' + os.path.join(TMP_FOLDER,TARBALL_FOLDER) + ' Exception: ' + str(sys.exc_info()[1])+'. Will try to continue, but it may fail.')
   
   logger.info(' creating bash script to run dq2-put ')
   # now write a script that setups up athena and DQ2 and uploads the dataset to the site
   out = ("""#!/usr/bin/env bash
# setup ASetup environment
export ATLAS_LOCAL_ROOT_BASE=""" + ATLAS_LOCAL_ROOT_BASE + """
# setup ASetup environment
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
# setup dq2 tools
echo ' setting up DQ2Client '
localSetupDQ2Client --skipConfirm
# upload dataset
echo ' running dq2-put to upload data set '
dq2-put -a -L """ + site + ' -s ' + os.path.join(TMP_FOLDER,TARBALL_FOLDER) + ' ' + outDS + '\n'
         )
   
   file = open(TMP_SHELL_SCRIPT,'w')
   file.write(out)
   file.close()
   file = None
   
   # run script
   logger.info(' running dq2-put bash script ')
   p = subprocess.Popen('sh ' + TMP_SHELL_SCRIPT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
   line = p.stdout.readline()
   while line != '':
      logger.info(line[0:-1])
      line = p.stdout.readline()
   p.wait()
   if p.returncode != 0:
      stdout,stderr = p.communicate()
      logger.error('ERROR Setting up athena and running prun command. returncode = ' + str(p.returncode) + '    stderr: \n' + stderr)
      return
   logger.debug(' finished creating dataset for alpgen output ' )



if __name__ == "__main__":
   sys.exit(main())


