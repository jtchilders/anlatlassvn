#!/usr/bin/env python
import os,sys,optparse,logging,subprocess
logger = logging.getLogger(__name__)

TMP_ATHENA_SCRIPT_FILENAME='tmp_athena_script.sh'

def main():
   logging.basicConfig(level=logging.INFO)
   parser = optparse.OptionParser(description='Submit athena condor job')
   parser.add_option('-g','--run-gen-trf',dest='run_gen_trf',help='Run Generate_trf instead of athena.',action='store_true',default=False)
   parser.add_option('-a','--run-atlasG4-trf',dest='run_atlasG4_trf',help='Run AtlasG4_trf instead of athena.',action='store_true',default=False)
   parser.add_option('-d','--run-digi-trf',dest='run_digi_trf',help='Run Digi_trf instead of athena.',action='store_true',default=False)
   parser.add_option('-r','--run-reco-trf',dest='run_reco_trf',help='Run Reco_trf instead of athena.',action='store_true',default=False)
   parser.add_option('-x','--args',dest='args',help='Arguments to pass to athena or trf.',default='')
   parser.add_option('-j','--asetup-args',dest='asetup_args',help='Arguments used to setup Athena environment',default='')
   parser.add_option('-k','--atlas-local',dest='atlas_local_root_base',help='Path to ATLAS Local setup, typically called ATLAS_LOCAL_ROOT_BASE.',default='/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase')
   options,args = parser.parse_args()

   trf_count = 0
   executable = 'athena'
   if options.run_gen_trf:
      executable = 'Generate_trf.py'
      trf_count += 1
   if options.run_atlasG4_trf: 
      executable = 'AtlasG4_trf.py'
      trf_count += 1
   if options.run_digi_trf:
      executable = 'Digi_trf.py'
      trf_count += 1
   if options.run_reco_trf:
      executable = 'Reco_trf.py'
      trf_count += 1
   if trf_count > 1:
      parser.error(' options -g -a -d -r are mutually exclusive, only one can be specified. ')
   
   return run_athena_in_shell( 
                        executable,
                        options.args,
                        options.asetup_args,
                        options.atlas_local_root_base,
                      )


def run_athena_in_shell(
                         executable             = 'athena',
                         args                   = '',
                         asetup_args            = '',
                         ATLAS_LOCAL_ROOT_BASE  = '/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase',
                       ):
   
   # create script to run athena
   script_content = ("""#!/usr/bin/env bash
# setup ASetup environment
export ATLAS_LOCAL_ROOT_BASE=""" + ATLAS_LOCAL_ROOT_BASE + """
# setup ASetup environment
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh 
# run asetup
ASETUP_ARGS='""" + asetup_args + """'
echo Running asetup with args: $ASETUP_ARGS
asetup $ASETUP_ARGS
RETVAL=$?
if [[ $RETVAL != 0 ]]; then
   echo Error running asetup
   exit -1
fi


# run athena/trf
EXE=""" + executable + """
ARG='""" + args + """'
echo Running $EXE $ARG
eval $EXE $ARG
RETVAL=$?
if [[ $RETVAL != 0 ]]; then
   echo ERROR running athena
   exit -2
fi
cd ..
"""
) 

   # write content to script file
   script_file = open(TMP_ATHENA_SCRIPT_FILENAME,'w')
   script_file.write(script_content)
   script_file.close()
   script_file = None
   
   p = subprocess.Popen('sh ' + TMP_ATHENA_SCRIPT_FILENAME,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
   line = p.stdout.readline()
   while line != '':
      logger.info(line[0:-1])
      line = p.stdout.readline()

   p.wait()

   if p.returncode != 0:
      logger.error(' Athena returned non-zero error code: stderr = ')
      for line in p.stderr.readlines():
         logger.error(' ' + line[0:-1])
      return -1

   return 0
   



if __name__ == '__main__':
   sys.exit(main())


