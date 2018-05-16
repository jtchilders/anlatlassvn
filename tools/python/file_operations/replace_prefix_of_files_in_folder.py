#!/usr/bin/env python

import sys,glob,os,shutil,logging,optparse
logger = logging.getLogger(__name__)

def main():
   logging.basicConfig(level=logging.INFO)
   parser = optparse.OptionParser(description="Globs files with the input prefix, copies them to files using the new prefix, and deletes the original files (optional).")
   parser.add_option('-i','--input-prefix',dest='input_prefix',help='Prefix of the files that will be changed.')
   parser.add_option('-o','--output-prefix',dest='output_prefix',help='The new prefix that will replace the input prefix')
   parser.add_option('-n','--no-delete',dest='delete',help='Do not delete the input files. Keep files with both prefixes.',action='store_false',default=True)
   parser.add_option('-f','--force-overwrite',dest='force_overwrite',help='If the destination file already exists, overwrite it.',action='store_true',default=False)
   opts,args = parser.parse_args()

   if opts.input_prefix is None or opts.output_prefix is None:
      parser.error(' must specify input and ouput prefix ')

   replace_prefix_of_files(opts.input_prefix,opts.output_prefix,opts.delete,opts.force_overwrite)

def replace_prefix_of_files(input_prefix,output_prefix,delete=True,force_overwrite=False):
   input_filenames = glob.glob(input_prefix + '*')
   for input_filename in input_filenames:
      postfix = input_filename[len(input_prefix):]
      output_filename = output_prefix + postfix

      logger.info('copying: %s -> %s' % (input_filename,output_filename))
      if os.path.exists(output_filename):
         if force_overwrite:
            logger.info('   destination file exists force_overwrite flag set so overwriting it. ')
            shutil.copyfile(input_filename,output_filename)
            if delete:
               logger.debug('    deleting input file ')
               os.remove(input_filename)
         else:
            logger.warning('   destination file exists, force_overwrite flag not set so skipping file.')
      else:
         shutil.copyfile(input_filename,output_filename)
         if delete:
            logger.debug('    deleting input file ')
            os.remove(input_filename)



if __name__ == '__main__':
   main()
