#!/usr/bin/env python
import os,sys,optparse,logging,subprocess,glob,time
logger = logging.getLogger(__name__)

TMP_FOLDER='/tmp/gen_trf_dagman/'
DAGMAN_FILENAME ='alpgen_gentrf.dag'
BASEJOB_FILENAME='alpgen_gentrf.condor'
BASEJOB_CONTENT = '''Universe                 = vanilla
Arguments                = -e $(ecmEnergy) -o $(outputEVNTFile) -i $(inputGeneratorFileArg) -j $(evgenJobOpts) -c $(jobConfig) -r $(runNumber)
Output                   = job.$(jobID).stdout.txt
Error                    = job.$(jobID).stderr.txt
Log                      = job.$(jobID).condor.txt
should_transfer_files    = YES
when_to_transfer_output  = ON_EXIT
Executable               = /users/hpcusers/svn/tools/python/athena/run_gen_trf_with_alpgen_input.py
transfer_input_files     = $(inputGeneratorFilePath)
transfer_output_files    = $(outputEVNTFile)
environment              = PYTHONPATH=/users/hpcusers/svn/tools/python
Queue'''


def main():
   logging.basicConfig(level=logging.INFO)
   parser = optparse.OptionParser(description='Create a Condor DAGMAN job to run many subjobs of Generate_trf on alpgen.')
   parser.add_option('-e','--ecmEnergy',dest='ecmEnergy',help='Center-of-mass energy.',type=int)
   parser.add_option('-o','--outputEVNTFile',dest='outputEVNTFile_base',help='Output EVNT file name base.')
   parser.add_option('-i','--input-glob',dest='input_glob',help='A command line GLOB with which the input tarballs will be located, for example "/path/to/tarballs/*.tgz".')
   parser.add_option('-j','--evgenJobOpts',dest='evgenJobOpts',help='Event Generation Job Options Tarball name.')
   parser.add_option('-c','--jobConfig',dest='jobConfig',help='Job Config (job options file).')
   parser.add_option('-r','--runNumber',dest='runNumber',help='Run Number.',type=int)
   parser.add_option('-f','--folder',dest='foldername',help='Folder name where dagman files will be written')
   options,args = parser.parse_args()

   if options.ecmEnergy is None:
      parser.error('must specify center of mass energy.')
   if options.input_glob is None:
      parser.error('must specify output EVNT file name base.')
   if options.jobConfig is None:
      parser.error('must specify job config file.')
   if options.runNumber is None:
      parser.error('must specify run number.')
   if options.evgenJobOpts is None:
      parser.error('must specify event generation job options.')
   if options.outputEVNTFile_base is None:
      parser.error('must specify event generation job options.')

   folder_name = None
   if options.foldername is not None:
      folder_name = options.foldername
   else:
      # create tmp folder for job files
      folder_name = os.path.join(TMP_FOLDER,str(int(time.time()*1000000)))
      logger.info(' using working path: ' + working_path)
   if not os.path.exists(folder_name):
      os.makedirs(folder_name)

   return create_alpgen_gen_trf_dagman(
                                options.input_glob,
                                options.ecmEnergy,
                                options.jobConfig,
                                options.runNumber,
                                options.evgenJobOpts,
                                options.outputEVNTFile_base,
                                folder_name,
                               )

def create_alpgen_gen_trf_dagman(input_glob,
                                 ecmEnergy,
                                 jobConfig,
                                 runNumber,
                                 evgenJobOpts,
                                 outputEVNTFile_base,
                                 working_path,
                                ):
   
   
   # create the base job file
   base_job = open(os.path.join(working_path,BASEJOB_FILENAME),'w')
   base_job.write(BASEJOB_CONTENT)
   base_job.close()
   base_job = None

   # create the dagman
   input_tarballs = glob.glob(input_glob)
   if len(input_tarballs) == 0:
      logger.error('no input files for glob: ' + str(input_glob))
      return -1
   
   dag_job = open(os.path.join(working_path,DAGMAN_FILENAME),'w')
   counter = 0
   for tarball in input_tarballs:
      str_counter = '%08i' % counter
      # define job
      dag_job.write('Job  %s %s\n' % (str_counter,BASEJOB_FILENAME))
      dag_job.write('VARS %s ecmEnergy="%i"\n' % (str_counter,ecmEnergy))
      dag_job.write('VARS %s inputGeneratorFileArg="%s"\n' % (str_counter,os.path.basename(tarball)))
      dag_job.write('VARS %s inputGeneratorFilePath="%s"\n' % (str_counter,tarball))
      dag_job.write('VARS %s jobConfig="%s"\n' % (str_counter,jobConfig))
      dag_job.write('VARS %s runNumber="%i"\n' % (str_counter,runNumber))
      dag_job.write('VARS %s evgenJobOpts="%s"\n' % (str_counter,evgenJobOpts))
      dag_job.write('VARS %s jobID="%s"\n' % (str_counter,str_counter))
      dag_job.write('VARS %s outputEVNTFile="%s.%s.pool"\n' % (str_counter,outputEVNTFile_base,str_counter))
      dag_job.write('\n')
      counter += 1
   dag_job.close()

   return 0


if __name__ == '__main__':
   sys.exit(main())

