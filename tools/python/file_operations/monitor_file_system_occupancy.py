#!/usr/bin/env python
import os,sys,optparse,logging,subprocess
sys.path.append('/users/hpcusers/svn/tools/python/')
import Mail
logger = logging.getLogger(__name__)

def main():
   logging.basicConfig(level=logging.INFO)

   parser = optparse.OptionParser(description='')
   parser.add_option('-p','--path',dest='path',help='File path to check.')
   parser.add_option('-e','--email',dest='email',help='The e-mail address to send the alert.')
   parser.add_option('-f','--fraction',dest='fraction',help='The occupancy fraction at which to send an alert.',type='float',default=0.9)
   parser.add_option('-a','--tmp-file',dest='tmp_file',help='The name of the temporary file where the current disk usage level will be stored and used to check if the disk usage has changed. That way you do not receive many warnings about an unchanged usage.')
   options,args = parser.parse_args()

   
   manditory_args = [
                     'path',
                     'email',
                     'fraction',
                  ]

   for man in manditory_args:
      if not options.__dict__[man]:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         return -1

   if options.tmp_file is None:
      options.tmp_file = os.path.join('/tmp/','monitor_usage_'+os.path.basename(options.path).replace('.','_')+'.log')
   logger.info('temporary file: ' + options.tmp_file)

   cmd = 'df ' + options.path
   p = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   stdout,stderr = p.communicate()
   fraction = 0.
   if p.returncode == 0:
      # first line of output is 'Filesystem      Size  Used Avail Use% Mounted on'
      # second line is          '-                11T  8.9T  1.5T  86% /grid/atlas/hpc'
      try:
         lines = stdout.split('\n')
         parts = lines[1].split()
         percent = parts[4]
         fraction = int(percent.replace('%',''))/100.
      except:
         logger.exception('Failed to parse "df" output: ' + stdout + '\n' + stderr + '\n')
         return -1
   else:
      logger.error('Failed to parse "df" output: ' + stdout + '\n' + stderr + '\n')
      return -1
   logger.info( 'fraction: ' + str(fraction))
   if fraction > options.fraction:
      # check if tmp file exists
      old_fraction = None
      if os.path.exists(options.tmp_file):
         # read in old fraction
         f = open(options.tmp_file)
         l = f.readlines()
         if len(l) >= 1:
            try:
               old_fraction = float(l[0].replace('\n',''))
            except:
               logger.exception('could not read old fraction.')
         f.close()

      if old_fraction is None or (old_fraction is not None and fraction > old_fraction):
         body = ' Path: ' + options.path + '\n'
         body += ' Partition Occupancy: ' + str(fraction) + '\n'
         Mail.send_mail('jchilders@anl.gov',options.email,'Monitoring ' + options.path,body)
         logger.info('email sent')

      # write fraction to tmp file
      f = open(options.tmp_file,'w')
      f.write(str(fraction))
      f.close()

   return 0

   


if __name__ == "__main__":
   sys.exit(main())