#!/usr/bin/env python
import os,sys,optparse,logging,subprocess
from mysql.mysql_wrapper import MySQL
logger = logging.getLogger(__name__)

DEFAULT_TABLE_NAME='argo_core_argodbentry'
DEFAULT_SERVER_NAME='localhost'
DEFAULT_USERNAME='root'
DEFAULT_PASSWORD=''
DEFAULT_DB_NAME='argo_cluster_production'
DEFAULT_RUCIO_SCOPE='group.phys-gener'
DEFAULT_RUCIO_SITE='ANLASC_SCRATCHDISK'

DEFAULT_SCOPE='group.phys-gener'
DEFAULT_DS_PREFIX='group.phys-gener.sherpa020101.'

SHERPA_BUILD_SCRIPT='/users/hpcusers/svn/tools/python/sherpa/sherpa_integration_presubmit.sh'
ARGO_SHERPA_TARBALL_NAME='sherpa_integration_output.tar.gz'

WORKING_FOLDER='/tmp/sherpa_tarball'

TARBALL_POSTFIX = '._00004.tar.gz'


UNNEEDED_FILES= [
                  'Process/*/P2_*',
                  'Sherpa_References.tex',
                  'sherpa.sh_checksum',
                  'job_submit.conf',
                  'sherpa_codegen_output.tar.gz',
                  'SConstruct',
                  'makelibs',
                  'run.sh.*',
                  'condor_*',
                ]
NEEDED_FILES = [
                'integration.log',
                'codegen.log',
               ]

def main():
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='upload sherpa jobs from ARGO to the grid')
   parser.add_option('-g','--group-pattern',dest='group_pattern',help='A pattern to use when searching the database by "group_identifier" to get jobs and upload them to the grid.')
   parser.add_option('-p','--process-folder',dest='process_folder',help='location of Process folder with libraries to copy into tar ball.')
   options,args = parser.parse_args()

   
   manditory_args = [
                     'group_pattern',
                     'process_folder',
                  ]

   for man in manditory_args:
      if options.__dict__[man] is None:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)

   rucio_upload_integration(options.group_pattern,options.process_folder)
   
def rucio_upload_integration(group_pattern,
                             process_folder,
                             db_server_name        = DEFAULT_SERVER_NAME,
                             db_server_login       = DEFAULT_USERNAME,
                             db_server_password    = DEFAULT_PASSWORD,
                             db_name               = DEFAULT_DB_NAME,
                             db_table_name         = DEFAULT_TABLE_NAME,
                             rucio_scope           = DEFAULT_RUCIO_SCOPE,
                             rucio_site            = DEFAULT_RUCIO_SITE,
                            ):
   
   mysql = MySQL(db_server_name,db_server_login,db_server_password,db_name)
   
   logger.info('Going to Upload using group_pattern = ' + group_pattern )
   logger.info('Using Process folder: ' + process_folder)

   where_conditions = "group_identifier LIKE '" + group_pattern +  "' AND state_current='HISTORY'"

   try:
      entries,labels = mysql.select(where_conditions=where_conditions)
   except Exception,e:
      logger.erro('mysql select command failed:' + str(e))
      raise

   logger.info('retrieved ' + str(len(entries)) + ' from DB.')

   # move to working folder
   start_folder = os.getcwd()
   if not os.path.exists(WORKING_FOLDER):
      os.mkdir(WORKING_FOLDER)
   os.chdir(WORKING_FOLDER)

   logger.info(' working folder: ' + WORKING_FOLDER)
   
   n = 0
   n_total = str(len(entries))
   for entry in entries:
      n += 1
      logger.info(' uploading ' + str(n) + ' of ' + n_total + ' jobs.')
      id = entry[labels['id']]
      group_identifier = entry[labels['group_identifier']]
      #logger.info(' Entry id = ' + str(id) + '  group_identifier = ' + group_identifier + '  state_current = ' + entry[label['state_current']] )

      job_dir = entry[labels['output_url']]
      if '://' in job_dir:
         job_dir = job_dir.split('://')[-1] # remove protocol
         job_dir = job_dir[job_dir.find('/'):] # remove server name

      logger.info('  job directory = ' + job_dir)
      
      # copy job files to the temproary folder
      for file in NEEDED_FILES:
         cmd = 'cp -r ' + job_dir + '/' + file + ' ' + WORKING_FOLDER + '/'
         logger.info(cmd)
         try:
            run(cmd,True)
         except Exception,e:
            logger.error('[id=' + str(entry[labels['id']]) + '] error in cp command: ' + str(e))
            raise

      sherpa_tarball = ARGO_SHERPA_TARBALL_NAME
      if os.path.exists(os.path.join(job_dir,sherpa_tarball)):
         cmd = 'cp -r ' + os.path.join(job_dir,sherpa_tarball) + ' ' + WORKING_FOLDER + '/'
         logger.info(cmd)
         try:
            run(cmd,True)
         except Exception,e:
            logger.error('[id=' + str(entry[labels['id']]) + '] error in cp command: ' + str(e))
            raise

         # extract DB from tarball and remove tarball
         cmd = 'tar zxvf ' + sherpa_tarball + ' Results.db'
         logger.info(cmd)
         try:
            run(cmd,True)
         except Exception,e:
            logger.error('[id=' + str(entry[labels['id']]) + '] Error running tar: ' + str(e))
            raise
         cmd = 'rm -rf ' + sherpa_tarball
         logger.info(cmd)
         try:
            run(cmd,True)
         except Exception,e:
            logger.error('[id=' + str(entry[labels['id']]) + '] Error running rm: ' + str(e))
            raise
      elif os.path.exists(os.path.join(job_dir,'Results.db')):
         cmd = 'cp -r ' + os.path.join(job_dir,'Results.db') + ' ' + WORKING_FOLDER + '/'
         logger.info(cmd)
         try:
            run(cmd,True)
         except Exception,e:
            logger.error('[id=' + str(entry[labels['id']]) + '] error in cp command: ' + str(e))
            raise
      else:
         logger.error(' No Results.db file found!')
         raise Exception('No Results.db file found!')

      # copy Process folder over from previous run
      cmd = 'cp -r ' + process_folder + ' ' + WORKING_FOLDER + '/'
      logger.info(cmd)
      try:
         run(cmd,True)
      except Exception,e:
         logger.error('[id=' + str(entry[labels['id']]) + '] Error running cp: ' + str(e))
         raise

      # get output file name
      index = group_identifier.find('_i1')
      container_name = group_identifier
      if index >= 0:
         container_name = group_identifier[:index]
      logger.info('  container: ' + container_name)

      dataset_name = group_identifier

      tarball_name = container_name + TARBALL_POSTFIX

      # tar needed files
      cmd = 'tar zcf ' + tarball_name + ' Process Results.db integration.log codegen.log'
      logger.info(cmd)
      try:
         run(cmd,True)
      except Exception,e:
         logger.error('[id=' + str(entry[labels['id']]) + '] Error running tar: ' + str(e))
         raise
      
      # now rucio upload
      cmd = 'rucio upload --scope ' + rucio_scope + ' --rse ' + rucio_site + ' ' + rucio_scope + ':' + dataset_name + ' ' + tarball_name
      logger.info(cmd)
      try:
         run(cmd)
      except Exception,e:
         if 'already exists on RSE. Will not try to reupload' in str(e):
            logger.info(' tarball already uploaded. Trying to continue')
         else:
            logger.error('[id=' + str(entry[labels['id']]) + '] Exception while trying to upload tarball: ' + str(e))
            raise


      # rucio create container
      cmd = 'rucio add-container ' + rucio_scope + ':' + container_name
      logger.info(cmd)
      try:
         run(cmd)
      except Exception,e:
         if 'Data Identifier Already Exists' in str(e):
            pass
         else:
            logger.error('[id=' + str(entry[labels['id']]) + '] Error in rucio: ' + str(e))
            raise

      # rucio add dataset to container
      cmd = 'rucio attach ' + rucio_scope + ':' + container_name + ' ' + rucio_scope + ':' + dataset_name
      logger.info(cmd)
      try:
         run(cmd)
      except Exception,e:
         if 'Data identifier already added to the destination content' in str(e):
            pass
         else:
            raise

      os.system('rm -rf ' + WORKING_FOLDER + '/*')



   os.chdir(start_folder)
   os.system('rm -rf ' + WORKING_FOLDER )
   logger.info('Done.')

def run(cmd,shell=False):
   tmp_cmd = cmd.split()
   if shell:
      tmp_cmd = cmd
   p = subprocess.Popen(tmp_cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=shell)
   lines = ''
   for line in p.stdout:
      logger.info(line[:-1])
      lines += line
   stdout,stderr = p.communicate()
   if p.returncode != 0:
      raise Exception(' cmd: ' + cmd + ' failed with exit code: ' + str(p.returncode) + '\n' + lines)

   return lines


    
def get_container_name(runnumber,descr):
   container_name = DEFAULT_DS_PREFIX + runnumber + '.' + descr.replace('_BFilter_','_') + '_HPC.TXT.mc15_v1'
   return container_name
   

def get_column_labels(cursor):
   description = cursor.description
   label = {}
   for i in range(len(description)):
      label[description[i][0]] = i
   return label



if __name__ == "__main__":
   main()
