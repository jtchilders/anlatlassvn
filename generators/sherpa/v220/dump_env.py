#!/usr/bin/env python
import os,sys,optparse,logging,subprocess
logger = logging.getLogger(__name__)

def main():
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='')
   parser.add_option('-o','--output',dest='output',help='output file+path')
   options,args = parser.parse_args()

   
   manditory_args = [
                     'output',
                  ]

   for man in manditory_args:
      if options.__dict__[man] is None:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)
   
   # get environment
   p = subprocess.Popen('env',stdout=subprocess.PIPE)
   stdout,stderr = p.communicate()

   env = stdout

   # get sherpa version




if __name__ == "__main__":
   main()
