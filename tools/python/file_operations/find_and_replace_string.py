#!/usr/bin/env python
import os,sys,optparse,logging,subprocess,shutil
logger = logging.getLogger(__name__)

def main():
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='')
   parser.add_option('-i','--str-to-replace',dest='string_to_replace',help='grep for this string in files.')
   parser.add_option('-o','--replace-with-str',dest='replace_with_string',help='replace grep-ed string with this string.')
   parser.add_option('-f','--replace-in-file',dest='filename',help='Replace only in this file. Cannot be used with "replace-in-path" option.')
   parser.add_option('-p','--replace-in-path',dest='path',help='Replace in this path.Cannot be used with "replace-in-file" option')
   parser.add_option('-r','--recursive',dest='recursive',help='If "repalce-in-path" option given, replace in all subdirectories as well.',action='store_true',default=False)
   options,args = parser.parse_args()

   
   manditory_args = [
                     'string_to_replace',
                     'replace_with_string',
                  ]

   for man in manditory_args:
      if options.__dict__[man] is None:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)

   if ((options.filename is not None and options.path is not None) or 
       (options.filename is None and options.path is None)):
      logger.error('Must specfiy one and only one of "replace-in-file" and "replace-in-path".')
      parser.print_help()
      sys.exit(-1)

   if options.recursive and options.path is None:
      logger.error('"recursive" option is set, but "replace-in-path" option is not.')
      parser.print_help()
      sys.exit(-1)

   # grep for files:

   grep_cmd = 'grep -l -I'
   if options.recursive:
      grep_cmd += ' -r'
   grep_cmd += ' "' + options.string_to_replace + '" '
   if options.filename is not None:
      grep_cmd += ' ' + options.filename
   elif options.path is not None:
      grep_cmd += ' ' + options.path
      if grep_cmd[-1] != '/' and grep_cmd[-1] != ' ':
         grep_cmd += '/'

   logger.info('running ' + grep_cmd) 

   p = subprocess.Popen(grep_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
   stdout,stderr = p.communicate()
   if p.returncode != 0:
      logger.error('grep returned not zero returncode: ' + str(p.returncode) + ': stdout = \n' + stdout + '\n stderr = \n' + stderr)
      sys.exit(-2)

   filelist = stdout.split('\n')
   

   escaped_string_to_replace = options.string_to_replace.replace('\\','\\\\').replace('/','\/').replace('.','\.')
   escaped_replace_with_string = options.replace_with_string.replace('\\','\\\\').replace('/','\/').replace('.','\.')

   for filename in filelist:
      if os.path.exists(filename):
         sed_cmd = 'sed "s/' + escaped_string_to_replace + '/' + escaped_replace_with_string + '/g" ' + filename
         logger.info('running ' + sed_cmd)
         p = subprocess.Popen(sed_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
         stdout,stderr = p.communicate()
         if p.returncode != 0:
            logger.error('sed returned not zero returncode: ' + str(p.returncode) + ': stdout = \n' + stdout + '\n stderr = \n' + stderr)
            sys.exit(-2)
         f = open( filename + '.tmp','w')
         f.write(stdout)
         f.close()
         shutil.move(filename+'.tmp',filename)


if __name__ == "__main__":
   main()
