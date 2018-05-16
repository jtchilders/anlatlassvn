#!/usr/bin/env python
import os,sys,optparse,logging,subprocess,shlex,glob,math
logging.basicConfig(format='%(asctime)s %(levelname)s:%(name)s:%(message)s',level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
   parser = optparse.OptionParser(description='loop over pool files, run checkFile.py and count up all the events')
   parser.add_option('-i','--input-glob',dest='input_glob',help='GLOB string to use to get a list of the files to check')
   options,args = parser.parse_args()
   
   if options.input_glob is None:
      parser.error('must specify -i')

   logger.info( ' input glob: ' + options.input_glob )
   filenames = glob.glob(options.input_glob)

   event_count = 0
   event_count2 = 0
   file_count = 0
   print 'looping over ' + str(len(filenames)) + ' files'
   for filename in filenames:
      file_count += 1
      p = subprocess.Popen(['checkFile.py','--fast',filename],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      line = p.stdout.readline()
      while line != '':
         logger.info(line[0:-1])
         
         if line.find('Nbr Events:') >= 0:
            split = line[0:-1].split()
            nevt = int(split[2])
            event_count += nevt
            event_count2 += nevt*nevt
            break
         line = p.stdout.readline()

      p.wait()
      stdout,stderr = p.communicate()
      if p.returncode != 0:
         logger.error('Error running checkFile.py, returncode = ' + str(p.returncode) + '\n stderr \n' + stderr)
         continue
      
      logger.info( ' >>>>>  events counted: ' + str(event_count) + ' files counted: ' + str(file_count)  )

   logger.info('event total: ' + str(event_count))
   
   if file_count != 0:
      mean = event_count / file_count
      logger.info(' average events per file: ' + str(mean))
      sigma = math.sqrt( (1./file_count)*event_count2 - mean*mean )
      logger.info(' spread on average: ' + str(sigma))


if __name__ == '__main__':
   main()



