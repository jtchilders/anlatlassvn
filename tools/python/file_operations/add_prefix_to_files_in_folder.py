#!/usr/bin/env python

import sys,glob,os,shutil,logging,optparse
logger = logging.getLogger('__main__')

def main():
   parser = optparse.OptionParser(description=" adds  a prefix to all the files in a folder")
   parser.add_option('-f','--folder',dest='folder',help='folder where files are located')
   parser.add_option('-p','--prefix',dest='prefix',help='prefix to add to files')
   parser.add_option('-o','--output-path',dest='output_path',help='optional ouput path, otherwise files go back into input folder')
   opts,args = parser.parse_args()

   if opts.folder is None or opts.prefix is None:
      parser.error(' must specify all options ')

   add_prefix_to_files_in_folder(opts.folder,opts.prefix,opts.output_path)

def add_prefix_to_files_in_folder(folder,prefix,output_path=None):
   files = glob.glob(folder + '/*')
   for file in files:
      old = file
      new = os.path.join(folder,prefix+os.path.basename(file))
      if output_path is not None:
         new = os.path.join(output_path,prefix+os.path.basename(file))
      logger.info('copying: ' + ('%50s' % old) + ' -> ' + ('%50s' % new))
      shutil.copyfile(old,new)





if __name__ == '__main__':
   main()
