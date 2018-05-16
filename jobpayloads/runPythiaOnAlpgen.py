#!/usr/bin/env python

import logging
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')
logger = logging.getLogger(__name__)
import optparse,os,sys,time,shutil,subprocess

from CondorJob import CondorJob

PYTHIA_EXE              = '/users/hpcusers/svn/generators/pythia/v8180/usercode/alpToHepmc/trunk/runPythiaOnAlpgen'
PYTHIA_INSTALL          = '/users/hpcusers/svn/generators/pythia/v8180/bgq/trunk'
PYTHIA_HEPMC_FILEBASE   = 'pythia_output.'
PYTHIA_MERGE_HEPMC      = '/users/hpcusers/balsam/exe/merge_hepmc.py'

def main():
   # print command
   print str(sys.argv)
   # parse command line
   parser = optparse.OptionParser(description='submit condor jobs for running pythia on alpgen output files')
   parser.add_option('-i','--input-base',dest='input_base',help='base string for input files')
   parser.add_option('-n','--nfiles',dest='nfiles',help='number of files to process')

   options,args = parser.parse_args()
   
   for i in range(options.nfiles):
      num_str = ('%08d' % i)
      base = options.input_base + num_str
      output = PYTHIA_HEPMC_FILEBASE + num_str + '.hepmc'
      args  = ' -x ' + os.path.join(PYTHIA_INSTALL,'xmldoc')
      args += ' -o ' + output
      args += ' -a ' + base
      job = CondorJob(executable = PYTHIA_EXE,
                       transfer_input_files     = (base  + '.unw,'
                                                  +base + '_unw.par,'
                                                  +PYTHIA_MERGE_HEPMC
                                                  ),
                       transfer_output_files    = output,
                       postcmd                  = os.path.basename(PYTHIA_MERGE_HEPMC),
                       output                   = 'pythia_' + num_str + '.stdout',
                       error                    = 'pythia_' + num_str + '.stderr',
                       log                      = 'pythia_' + num_str + '.condor_log',
                       arguments                = args,
                      )
      job.submit()

   
   
if __name__ == '__main__':
   main()


