#!/usr/bin/env python
import os,sys,optparse,logging,shutil,glob,subprocess,stat
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

JOB_OPTS = os.path.join(os.path.dirname(os.path.realpath(__file__)),'ConvertHepMcToPool_jobOptions.py')
TMP_JOB_OPTS = 'tmp_jo.py'
INSERT_INPUT_FILE_HERE = 'INSERT_INPUT_FILE_HERE'
INSERT_OUTPUT_FILE_HERE = 'INSERT_OUTPUT_FILE_HERE'

# need to setup Athena
ATLAS_LOCAL_ROOT_BASE = '/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase'
ASETUP_ARGS = ' --testarea=/users/jchilders/athenaTestarea/17.7.3.6  17.7.3.6,AtlasProduction,slc5'
ATHENA_SCRIPT_FILENAME = 'athena.sh'
ATHENA_SCRIPT_CONTENT = """#!/usr/bin/env bash
# setup ASetup environment
export ATLAS_LOCAL_ROOT_BASE=""" + ATLAS_LOCAL_ROOT_BASE + """
# setup ASetup environment
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
echo ' setting up asetup '
asetup """ + ASETUP_ARGS + """
echo ' running athena '
athena -c "$1" """ + JOB_OPTS + "\n"

DAGMAN_FILENAME ='convertHepMcToPool.dag'
DAGMAN_CONTENT = '''Job %08i %s
VARS %08i jobID="%08i"
VARS %08i inputHepmc="%s"
VARS %08i outputEVNTFile="%s"
VARS %08i athena_args="%s"

'''  
BASEJOB_FILENAME='convertHepMcToPool.condor'
BASEJOB_CONTENT = '''Universe                 = vanilla
Arguments                = sh ''' + ATHENA_SCRIPT_FILENAME + ''' $(athena_args)
Output                   = job.$(jobID).stdout.txt
Error                    = job.$(jobID).stderr.txt
Log                      = job.$(jobID).condor.txt
should_transfer_files    = YES
when_to_transfer_output  = ON_EXIT
Executable               = /usr/bin/env
transfer_input_files     = $(inputHepmc),''' + ATHENA_SCRIPT_FILENAME + '''
transfer_output_files    = $(outputEVNTFile)
environment              = PYTHONPATH=/users/hpcusers/svn/tools/python
Queue'''

def main():

   parser = optparse.OptionParser(description='')
   parser.add_option('-g','--input-glob',dest='input_glob',help='input glob pattern for the files to convert. Must use quotation around pattern. For Example: "/path/to/files/*.hepmc"')
   parser.add_option('-o','--output-path',dest='output_path',help='output destination folder. If this is not specified, the output pool files will be written to the input files location.')
   parser.add_option('-r','--run-number',dest='run_number',help='The run number to set in the meta-data of the EVNT file.')
   parser.add_option('-p','--out-filename',dest='out_filename',help='Output filename base.',default='datasetfile')
   parser.add_option('-d','--dump-dagman',dest='dump_dagman',action='store_true',
      help='Instead of running locally, dump a dagman job and exit.',default = False)
   options,args = parser.parse_args()

   if options.input_glob is None:
      parser.error('Must specify -g')
   if options.run_number is None:
      parser.error('Must specify -r')
   
   # get file list
   filelist = glob.glob(options.input_glob)
   logger.info(' using glob pattern: ' + options.input_glob)
   logger.info(' processing ' + str(len(filelist)) + ' files ')
   if options.output_path is not None:
      logger.info( ' output destination: ' + options.output_path )

   ConvertManyHepMcFilesToPool(filelist,options.output_path,options.run_number,options.out_filename,options.dump_dagman)

def ConvertManyHepMcFilesToPool(filelist,output_path = None,run_number = None,out_filename = None,dump_dagman = False):
   global BASEJOB_FILENAME,DAGMAN_FILENAME,ATHENA_SCRIPT_FILENAME
   if output_path is not None:
      if not os.path.exists(output_path):
         logger.error(' output path does not exist.')
         return -1
      
      DAGMAN_FILENAME  = os.path.join(output_path,DAGMAN_FILENAME)
      ATHENA_SCRIPT_FILENAME = os.path.join(output_path,ATHENA_SCRIPT_FILENAME)


   # create athena job options
   athenaSH = open(ATHENA_SCRIPT_FILENAME,'w')
   athenaSH.write(ATHENA_SCRIPT_CONTENT)
   athenaSH.close()


   # create condor base job
   if dump_dagman:
      base = open(os.path.join(output_path,BASEJOB_FILENAME),'w')
      base.write(BASEJOB_CONTENT)
      base.close()

      # delete dagman if it exists
      if os.path.exists(DAGMAN_FILENAME):
         os.remove(DAGMAN_FILENAME)

   # loop over each file and run athena
   i = 0
   for file in filelist:
      #logger.info(' processing file: ' + file)
      
      # create output filename
      out_file = file + '.pool'
      if out_filename is not None:
         out_file = out_filename + ('%08i' % i) + '.pool'

      # create athena run script

      # first create the variable list to pass to the job options file
      athena_args = 'inputFilename=\\\\\\\"' + os.path.basename(file) + '\\\\\\\";outputFilename=\\\\\\\"' + out_file + '\\\\\\\"'
      if run_number is not None:
         athena_args += ';runNumber='+str(run_number)
      
      if dump_dagman:
         
         # append job to dagman
         dagman = open(DAGMAN_FILENAME,'a')
         dagman.write(DAGMAN_CONTENT % (i,BASEJOB_FILENAME,i,i,i,file,i,out_file,i,athena_args))
         dagman.close()



      else:
         p = subprocess.Popen('sh ' + ATHENA_SCRIPT_FILENAME + ' "' + athena_args + '"',stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=subprocess.PIPE,shell=True)
         line = p.stdout.readline()
         while line != '':
            logger.info(line[0:-1])
            line = p.stdout.readline()
         p.wait()
         if p.returncode != 0:
            stdout,stderr = p.communicate()
            logger.error(' input file ' + file + ' produced errors during running: stderr = \n' + stderr )

      i += 1


if __name__ == "__main__":
   main()
