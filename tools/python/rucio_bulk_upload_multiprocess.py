#!/usr/bin/env python
import os,sys,optparse,logging,multiprocessing,glob,subprocess
import rucio_bulk_upload
import MySQLdb
logger = logging.getLogger(__name__)

USERNAME='usatlas2'

class OneTask:
   def __init__(self,
                container_name,
                dataset_postfix,
               ):
      self.container_name        = container_name
      self.dataset_name          = self.container_name + dataset_postfix
      self.glob                  = self.container_name + '/' + self.dataset_name + '/' + self.container_name + '._*.tar.gz'



tasks = [ 
  
   OneTask('group.phys-gener.alpgen214.365633.AlpgenPythia_P2012_ktfac0p5_ZeebbNp3Incl_HPC.TXT.mc15_v1','_i11'),
   OneTask('group.phys-gener.alpgen214.365643.AlpgenPythia_P2012_ktfac0p5_ZmumubbNp3Incl_HPC.TXT.mc15_v1','_i11'),
   OneTask('group.phys-gener.alpgen214.365433.AlpgenPythia_P2012_ktfac2_ZeebbNp3Incl_HPC.TXT.mc15_v1','_i11'),
   OneTask('group.phys-gener.alpgen214.365443.AlpgenPythia_P2012_ktfac2_ZmumubbNp3Incl_HPC.TXT.mc15_v1','_i11'),
   OneTask('group.phys-gener.alpgen214.365233.AlpgenPythia_P2012_qfac0p5_ZeebbNp3Incl_HPC.TXT.mc15_v1','_i11'),
   OneTask('group.phys-gener.alpgen214.365243.AlpgenPythia_P2012_qfac0p5_ZmumubbNp3Incl_HPC.TXT.mc15_v1','_i11'),
   OneTask('group.phys-gener.alpgen214.365033.AlpgenPythia_P2012_qfac2_ZeebbNp3Incl_HPC.TXT.mc15_v1','_i11'),
   OneTask('group.phys-gener.alpgen214.365043.AlpgenPythia_P2012_qfac2_ZmumubbNp3Incl_HPC.TXT.mc15_v1','_i11')


]



def run_upload(task):
   
   logger.info(' Running task with: ')
   logger.info('    glob:             ' + str(task.glob))
   logger.info('    dataset_name:     ' + str(task.dataset_name))
   logger.info('    container_name:   ' + str(task.container_name))
   
   # get file list

   file_list = glob.glob(task.glob)
   logger.info('    number of files:  ' + str(len(file_list)))
   if len(file_list) > 0:
      
      rucio_bulk_upload.rucio_bulk_upload(
            file_list,
            task.dataset_name,
            task.container_name,
            username = USERNAME,
          )
   else:
      logger.info(' no files for ' + task.container_name)

   cmd = 'rucio list-files ' + task.container_name
   p = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
   stdout,stderr = p.communicate()

   start = stdout.rfind('Total files : ')
   end = stdout.find('\n',start)
   if start > 0 and end > 0:
      logger.info(' container: ' + task.container_name + ' has ' + stdout[start+13:end] + ' files.')
   
   logger.info(' done with: ' + str(task.container_name))


def main():
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(process)d %(levelname)s:%(name)s:%(message)s')
   logger.info('Beginning pool with ' + str(len(tasks)) + ' tasks.')
   pool = multiprocessing.Pool(5)   
   pool.map(run_upload,tasks)
   logger.info('All processes is Done')


if __name__ == "__main__":
   main()
