#!/usr/bin/env python
import os,sys,optparse,logging,glob,subprocess,time,math
logger = logging.getLogger(__name__)

def main():
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='')
   parser.add_option('-g','--glob',dest='glob',help='Glob to use to get list of files to upload to dataset. Using quotation marks.')
   parser.add_option('-c','--container-name',dest='container_name',help='Name of container to which to add dataset.')
   parser.add_option('-d','--dataset-postfix',dest='dataset_postfix',help='Postfix to add to container name to create dataset name.')
   parser.add_option('-s','--scope',dest='scope',help='Rucio --scope option',default='group.phys-gener')
   parser.add_option('-r','--rse',dest='resource',help='Rucio --rse option',default='ANLASC_SCRATCHDISK')
   parser.add_option('-n','--num-retries',dest='num_retries',help='Number of retries when uploading fails before calling it quits.',default=10,type='int')
   parser.add_option('-l','--log-file',dest='log_filename',help='log file where progress is written',default='rucio_upload.log')
   options,args = parser.parse_args()

   
   manditory_args = [
                     'glob',
                     'container_name',
                     'dataset_postfix',
                     'scope',
                     'resource',
                     'num_retries',
                     'log_filename',
                  ]

   for man in manditory_args:
      if options.__dict__[man] is None:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)
   dataset_name = options.container_name + options.dataset_postfix
   logger.info('glob:            ' + options.glob)
   logger.info('container name:  ' + options.container_name)
   logger.info('dataset name:    ' + dataset_name)
   logger.info('scope:           ' + options.scope)
   logger.info('resource:        ' + options.resource)
   logger.info('num_retries:     ' + str(options.num_retries))
   logger.info('log_filename:    ' + options.log_filename)

   file_list = glob.glob(options.glob)

   logger.info(' glob returned ' + str(len(file_list)) + ' files ')

   duration_sum = 0.
   duration_sum2 = 0.
   duration_count = 0
   
   for file in file_list:
      cmd = 'rucio upload --scope ' + options.scope + ' --rse ' + options.resource + ' ' + options.scope + ':' + dataset_name + ' ' + file
      logger.info(cmd)
      start_time = time.time()
      for retry in range(options.num_retries):
         p = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
         for line in p.stderr:
            logger.info(line[:-1])
         p.wait()
         if int(p.returncode) != 0:
            stdout,stderr = p.communicate()
            logger.error(' rucio upload command return code non-zero = ' + str(p.returncode) + '\n stdout: \n' + stdout + '\n stderr: \n' + stderr)
            logger.error(' retrying ' + str(retry) + ' of ' + str(options.num_retries))
         else:
            break
      duration = time.time() - start_time
      duration_sum += duration
      duration_sum2 += duration*duration
      duration_count += 1
      ave_dur = duration_sum / duration_count
      logger.info(str(duration_count) + ' ' + str(ave_dur) + ' ' + str(duration_sum2))
      sigma_dur = math.sqrt((1./duration_count)*duration_sum2 - ave_dur*ave_dur)

      estimated_time_left = (len(file_list) - duration_count)*ave_dur
      hour = int(estimated_time_left/60./60.)
      min  = int((estimated_time_left - hour*60.*60.)/60.)
      sec  = int(estimated_time_left - hour*60.*60. - min*60.)

      if (duration_count % 10) == 0: 
         logger.info('average duration = ' + str(ave_dur))
         logger.info('duration sigma =   ' + str(sigma_dur))
         logger.info('uploaded files ' + str(file) + ' of ' + str(len(file_list)))
         logger.info('estimated time left: %02i:%02i:%02i' % (hour,min,sec))
       
   
   cmd = 'rucio add-container '+ options.scope + ':' + options.container_name
   logger.info(cmd)
   p = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   for line in p.stderr:
      logger.info(line[:-1])
   p.wait()
   if int(p.returncode) != 0:
      stdout,stderr = p.communicate()
      logger.error(' rucio add-container command return code non-zero = ' + str(p.returncode) + '\n stdout: \n' + stdout + '\n stderr: \n' + stderr)
      sys.exit(-1)

   cmd = 'rucio attach '+ options.scope + ':' + options.container_name + ' ' + dataset_name
   logger.info(cmd)
   p = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   for line in p.stderr:
      logger.info(line[:-1])
   p.wait()
   if int(p.returncode) != 0:
      stdout,stderr = p.communicate()
      logger.error(' rucio attach command return code non-zero = ' + str(p.returncode) + '\n stdout: \n' + stdout + '\n stderr: \n' + stderr)
      sys.exit(-1)

   logger.info(' finished uploading ' + str(len(file_list)) + ' files ')
         


if __name__ == "__main__":
   main()
