#!/usr/bin/env python
import os,sys,optparse,logging,subprocess
from mysql.mysql_wrapper import MySQL
logger = logging.getLogger(__name__)


def main():
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='upload sherpa jobs from ARGO to the grid')
   parser.add_option('-p','--process-folder',dest='process_folder',help='location of Process folder with libraries to copy into tar ball.')
   options,args = parser.parse_args()

   
   manditory_args = [
                     'process_folder',
                  ]

   for man in manditory_args:
      if options.__dict__[man] is None:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)




if __name__ == "__main__":
   main()
