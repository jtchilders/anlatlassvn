#!/usr/bin/env python
import os,sys,optparse,logging
logger = logging.getLogger(__name__)
sys.path.append('/users/hpcusers/svn/tools/python/')
from CondorJob import CondorJob

ATHENA_SCRIPT_FILENAME='tmp_athena_script.sh'
CONDOR_JOB_FILENAME='condor_submit.txt'
ATLASLocalRootBase='/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase'

def main():
   parser = optparse.OptionParser(description='Submit athena condor job')
   parser.add_option('-g','--run-gen-trf',dest='run_gen_trf',help='Run Generate_trf instead of athena.',action='store_true',default=False)
   parser.add_option('-a','--run-atlasG4-trf',dest='run_atlasG4_trf',help='Run AtlasG4_trf instead of athena.',action='store_true',default=False)
   parser.add_option('-d','--run-digi-trf',dest='run_digi_trf',help='Run Digi_trf instead of athena.',action='store_true',default=False)
   parser.add_option('-r','--run-reco-trf',dest='run_reco_trf',help='Run Reco_trf instead of athena.',action='store_true',default=False)
   parser.add_option('-x','--args',dest='args',help='Arguments to pass to athena or trf.',default='')
   parser.add_option('-j','--asetup-args',dest='asetup_args',help='Arguments used to setup Athena environment',default='')
   parser.add_option('-k','--atlas-local',dest='atlas_local_root_base',help='Path to ATLAS Local setup, typically called ATLAS_LOCAL_ROOT_BASE.',default=ATLASLocalRootBase)
   parser.add_option('-i','--input-files',dest='input_files',help='Comma separated list of files to transfer into Condor job',default='')
   parser.add_option('-o','--output-files',dest='output_files',help='Comma separated list of files to transfer out of Conodr job',default='')
   parser.add_option('--script-name',dest='script_name',help='The name of the temporary bash script that is run by the condor job.[default="' + ATHENA_SCRIPT_FILENAME + '"]',default=ATHENA_SCRIPT_FILENAME)
   parser.add_option('--condor-job-filename',dest='condor_job_filename',help='The name of the condor job file produced that is submitted to condor_submit .[default="' + CONDOR_JOB_FILENAME + '"]',default=CONDOR_JOB_FILENAME)
   options,args = parser.parse_args()

   trf_count = 0
   executable = 'athena'
   if options.run_gen_trf:
      executable = 'Generate_tf.py'
      trf_count += 1
   if options.run_atlasG4_trf: 
      executable = 'AtlasG4_tf.py'
      trf_count += 1
   if options.run_digi_trf:
      executable = 'Digi_tf.py'
      trf_count += 1
   if options.run_reco_trf:
      executable = 'Reco_tf.py'
      trf_count += 1
   if trf_count > 1:
      parser.error(' options -g -a -d -r are mutually exclusive, only one can be specified. ')
   
   athena_condor_job(executable,
                     options.args,
                     options.asetup_args,
                     options.atlas_local_root_base,
                     options.input_files,
                     options.output_files,
                     options.tmp_script_name,
                     options.condor_job_filename,
                    )


def athena_condor_job(
                      executable             = 'athena',
                      args                   = '',
                      asetup_args            = '',
                      ATLAS_LOCAL_ROOT_BASE  = '/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase',
                      input_files            = '',
                      output_files           = '',
                      tmp_script_name        = ATHENA_SCRIPT_FILENAME,
                      condor_job_filename    = CONDOR_JOB_FILENAME,
                      output                 = '$(cluster).stdout',
                      error                  = '$(cluster).stderr',
                      log                    = '$(cluster).condorlog',
                      output_file_renames    = [],
                      no_submit              = False,
                      output_remaps          = None,
                     ):
   
   # create script to run athena
   script_content = ("""#!/usr/bin/env bash
# setup ASetup environment
export ATLAS_LOCAL_ROOT_BASE=""" + ATLAS_LOCAL_ROOT_BASE + """
# setup ASetup environment
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh 
# run asetup
asetup """ + asetup_args + """
# run athena
""" + executable + ' ' + args + '\n') 

   # write content to script file
   script_file = open(tmp_script_name,'w')
   script_file.write(script_content)
   script_file.close()
   script_file = None

   if len(input_files) > 0:
      input_files += ',' + tmp_script_name
   else:
      input_files = tmp_script_name

   job = CondorJob(
                   executable = '/usr/bin/env',
                   arguments = 'sh ' + tmp_script_name,
                   transfer_input_files      = input_files,
                   transfer_output_files     = output_files,
                   transfer_output_remaps    = output_remaps,
                   output = output,
                   error = error,
                   log = log,
                  )

   job.write(condor_job_filename)
   if no_submit:
      print str(job)
   else:
      job.submit()
                   
   



if __name__ == '__main__':
   sys.exit(main())


