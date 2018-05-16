#!/usr/bin/env python
import os,sys,logging,optparse,shutil
from athena.run_athena_in_shell import run_athena_in_shell

ATHENA_INSTALLAREA='/users/hpcusers/athenaInstallArea'
GENERATE_TRF = 'Generate_trf.py'
ASETUP_ARGS  = '--testarea='+ATHENA_INSTALLAREA+' AtlasProduction,17.7.3.6,gcc46'
TMP_WORK_DIR = 'tmp_get_alpgen_input_card'
INPUT_FILENAME_IMODE1 = 'input_card.mode_1.dat'
INPUT_FILENAME_IMODE2 = 'input_card.mode_2.dat'

def main():
   logging.basicConfig(level=logging.INFO)
   parser = optparse.OptionParser(description='Generate Alpgen Input Card from Athena inputs.')
   parser.add_option('-e','--ecmEnergy',dest='ecmEnergy',help='Center of Mass energy in GeV for run. [default=13000GeV]',type=int,default=13000)
   parser.add_option('-r','--runNumber',dest='runNumber',help='Run number.',type=int)
   parser.add_option('-o','--jobConfig',dest='jobConfig',help='Job Options file.')
   parser.add_option('-j','--evgenJobOpts',dest='evgenJobOpts',help='Event job options tarball which is retrieved from a webserver. They look like this "MC12JobOptions-00-00-00_v1.tar.gz".')
   parser.add_option('-k','--keep-files',dest='keep_files',help='Keep athena tmp files.',action='store_true',default=False)
   options,args = parser.parse_args()

   if options.runNumber is None:
      parser.error('must specify run number')
   if options.jobConfig is None:
      parser.error('must specify job config')
   if options.evgenJobOpts is None:
      parser.error('must specify evgen job options')

   get_alpgen_input_card(
                         ecmEnergy           = options.ecmEnergy,
                         runNumber           = options.runNumber,
                         jobConfig           = options.jobConfig,
                         evgenJobOpts        = options.evgenJobOpts,
                         keep_files          = options.keep_files,
                        )

   return 0


def get_alpgen_input_card(ecmEnergy          = 13000, #GeV
                          runNumber          = 123456,
                          jobConfig          = 'some_joboptions.py',
                          evgenJobOpts       = 'MC12JobOpts-00-11-59_v5.tar.gz',
                          keep_files         = False,
                         ):
   
   # Move to temporary work directory:
   CWD = os.getcwd()
   os.mkdir(TMP_WORK_DIR)
   os.chdir(TMP_WORK_DIR)

   # build generate TRF args from inputs
   # example: ecmEnergy=13000 randomSeed=122323 outputEVNTFile=myoutput.pool runNumber=181964 firstEvent=1 
   #          jobConfig=~/athenaTestarea/17.2.11.12/ap1/MC12.181964.AlpgenPythia_Auto_P2011C_Dil_10_6_ZtautauNp4.py 
   #          inputGeneratorFile=alpgen.000001.TXT.v1.tgz evgenJobOpts=MC12JobOpts-00-11-59_v5.tar.gz 
   #          --athenaopts=""--config-only=test.txt""'   
   dummy_outputEVNTFile = 'myoutput.pool'
   dummy_randomSeed = 131345
   dummy_firstEvent  = 1
   dummy_athenaopts  = '--config-only=test.txt'

   ### Create the dummy inputGeneratorFile
   dummy_inputGeneratorFile = 'alpgen.00001.TXT.v1.tar.gz'
   os.system('touch alpgen.00001.TXT.v1._00001.dat')
   os.system('touch alpgen.00001.TXT.v1._00001.events')
   os.system('tar zcf ' + dummy_inputGeneratorFile + ' alpgen.00001.TXT.v1._00001.dat alpgen.00001.TXT.v1._00001.events')
   os.remove('alpgen.00001.TXT.v1._00001.dat')
   os.remove('alpgen.00001.TXT.v1._00001.events')

   
   ARGS = ''
   ARGS += ' ecmEnergy='+str(ecmEnergy)
   ARGS += ' randomSeed='+str(dummy_randomSeed)
   ARGS += ' outputEVNTFile='+dummy_outputEVNTFile
   ARGS += ' runNumber='+str(runNumber)
   ARGS += ' firstEvent='+str(dummy_firstEvent)
   ARGS += ' jobConfig='+jobConfig
#   ARGS += ' inputGeneratorFile='+inputGeneratorFile
   ARGS += ' evgenJobOpts='+evgenJobOpts
   ARGS += ' --athenaopts="'+dummy_athenaopts+'"'
   

   run_athena_in_shell(GENERATE_TRF,ARGS,ASETUP_ARGS)

   shutil.move(INPUT_FILENAME_IMODE1,CWD)
   shutil.move(INPUT_FILENAME_IMODE2,CWD)
   os.chdir(CWD)
   if not keep_files:
      shutil.rmtree(TMP_WORK_DIR)


if __name__ == '__main__':
   sys.exit(main())


