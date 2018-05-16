#!/usr/bin/env python
import os,sys,optparse,logging,glob
import submit_sherpa
logger = logging.getLogger(__name__)

def main():
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='')
   parser.add_option('-g','--input-glob',dest='input_glob',help='The input pattern for finding job options. Use quotation marks.')
   options,args = parser.parse_args()

   
   manditory_args = [
                     'input_glob',
                  ]

   for man in manditory_args:
      if options.__dict__[man] is None:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)
   
   job_option_files = sorted(glob.glob(options.input_glob))

   if len(job_option_files) == 0:
      logger.error('Glob pattern failed to retrieve anything.')

   for job_option_file in job_option_files:
      logger.info(' Running over ' + job_option_file )

      # ./submit_sherpa.py  -b argo_cluster -e edison -f 3 -g 48 -i 600 -p sherpa.MC15.777022.Sherpa_CT10_Zee_Pt140_280_BFilter_qsf025.edison --python-input-card /grid/atlas/hpc/data/sherpa/inputcards/Zee/joboptions/MC15.777023.Sherpa_CT10_Zee_Pt140_280_BFilter_qsf4.py
      submit_sherpa.submit_sherpa(None,'argo_cluster',
                                  None,'edison',3,48,600,
                                  None,None,None,None,None,None,
                                  None,True,
                                  'sherpa.' + os.path.basename(job_option_file).replace('.py','') + '.edison',
                                  python_input_card = job_option_file
                                 )

   logger.info(' Done. ')
                               
                               

if __name__ == "__main__":
   main()
