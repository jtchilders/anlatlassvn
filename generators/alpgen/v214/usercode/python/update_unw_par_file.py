#!/usr/bin/env python
import optparse,os,shutil,logging,sys
logger = logging.getLogger(__name__)
from AlpgenUnwParFile import AlpgenUnwParFile

EXCLUSIVE_JET_MATCHING=1
INCLUSIVE_JET_MATCHING=0
def UpdateUnwParFile(input_filename,output_filename=None,jet_treatment=EXCLUSIVE_JET_MATCHING):
   
   if output_filename is None:
      output_filename = input_filename + '.tmp'

   parFile = AlpgenUnwParFile.read_file(input_filename)
   new_parameters_list = [
               AlpgenUnwParFile.AlpgenUnwParameter(501,20.,'min ETCLUS used for parton-jet matching (Normally ETCLUS = ptjmin + 5)'),
               AlpgenUnwParFile.AlpgenUnwParameter(502,0.7,'min RCLUS value for parton-jet matching'),
               AlpgenUnwParFile.AlpgenUnwParameter(503,6.0,'max ETACLUS value for parton-jet matching'),
               AlpgenUnwParFile.AlpgenUnwParameter(504,jet_treatment,'Jet Matching: 0 inclusive 1 exclusive'),
              ]
   parFile.parameters += new_parameters_list
   parFile.write_file(output_filename)
   try:
      shutil.copy(input_filename,input_filename+'.old')
      os.remove(input_filename)
      shutil.copy(output_filename,input_filename)
      os.remove(output_filename)
   except:
      logger.error(' exception caught: ' + str(sys.exc_info()[1]))
      raise


def main():
   parser = optparse.OptionParser(description='Add parameters to par file for Athena jobs')
   parser.add_option('-i','--input-file',dest='input_filename',help='input par file')
   parser.add_option('-o','--output-file',dest='output_filename',help='outut par file')
   parser.add_option('-j','--inclusive-jet',dest='inclusive_jet',help='If this flag is included, the jet treatment will be inclusive, otherwise it is exclusive.',action='store_true',default=False)
   options,args = parser.parse_args()


   if options.inclusive_jet:
      UpdateUnwParFile(options.input_filename,options.output_filename,INCLUSIVE_JET_MATCHING)
   else:
      UpdateUnwParFile(options.input_filename,options.output_filename,EXCLUSIVE_JET_MATCHING)

if __name__ == "__main__":
   main()
   


