#!/usr/bin/env python
import os,sys,optparse,logging,glob,subprocess
logger = logging.getLogger(__name__)

def main():
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='')
   parser.add_option('-a',dest='inputA',help='inputA, in case of duplicate this folder will be deleted')
   parser.add_option('-b',dest='inputB',help='inputB')
   options,args = parser.parse_args()

   
   manditory_args = [
                     'inputA',
                     'inputB',
                  ]

   for man in manditory_args:
      if options.__dict__[man] is None:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)
   
   
   folder_list_A = glob.glob(os.path.join(options.inputA + '*'))
   folder_list_B = glob.glob(os.path.join(options.inputB + '*'))
   
   common_folders = []

   for folderA in folder_list_A:
      common = False
      for folderB in folder_list_B:
         if os.path.basename(folderA) == os.path.basename(folderB):
            common_folders.append(os.path.basename(folderA))
            break
   print len(common_folders),'duplicate folders'
   for folder in common_folders:
      full_path = os.path.join(options.inputA,folder)
      print 'deleting',full_path
      os.system('rm -rf ' + full_path)
   

   '''
   duplicate_folders = []
   for folder in common_folders:
      cmd = 'diff -r -q ' + os.path.join(options.inputA,folder) + ' ' + os.path.join(options.inputB,folder)

      p = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
      stdout,stderr = p.communicate()
      print '"' + stdout + '"' 
      if p.returncode == 0 and len(stdout) == 0:
         print 'found duplicate: ',folder'''


if __name__ == "__main__":
   main()
