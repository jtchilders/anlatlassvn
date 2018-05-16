#!/usr/bin/env python

import os,sys,glob,optparse,logging
logger = logging.getLogger(__name__)

def main():
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='Build Alpgen Executables')
   parser.add_option('-g','--gen',dest='gen',help='Build the Alpgen Fortran 77 executable.',action='store_true',default=False)
   parser.add_option('-i','--gen90',dest='gen90',help='Build the Alpgen Fortran 90 executable. Tends to be faster execution time than Fortran 77 version.',action='store_true',default=False)
   parser.add_option('-m','--mpi',dest='mpi',help='Build Alpgen with MPI.',action='store_true',default=False)
   parser.add_option('-p','--cteq-only',dest='cteq_only',help='Build Alpgen with only CTEQ PDF. This reduces the executable size in memory. If running on Mira, this is required for reaching 64 ranks per node.',action='store_true',default=False)
   options,args = parser.parse_args()

   
   manditory_args = [
                     'gen',
                     'gen90',
                     'mpi',
                     'cteq_only',
                  ]

   for man in manditory_args:
      if options.__dict__[man] is None:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)

   if not options.gen and not options.gen90:
      logger.error('Must choose gen or gen90')
      parser.print_help()
      sys.exit(-1)

   # get a list of the make files
   file_list = sorted(glob.glob('*/Makefile'))
   # loop over them and compile the requested binaries
   for file in file_list:
      dir = file.split('/')[0]
      if options.mpi and options.cteq_only:
         if options.gen:
            os.system('make -C ' + dir + ' gen_mpi_nomrstpdfs')
         if options.gen90:
            os.system('make -C ' + dir + ' gen90_mpi_nomrstpdfs')
      elif options.mpi and not options.cteq_only:
         if options.gen:
            os.system('make -C ' + dir + ' gen_mpi')
         if options.gen90:
            os.system('make -C ' + dir + ' gen90_mpi')
      else:
         if options.gen:
            os.system('make -C ' + dir + ' gen')
         if options.gen90:
            os.system('make -C ' + dir + ' gen90')

if __name__ == "__main__":
   main()


