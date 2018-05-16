#!/usr/bin/env python
import os,sys,optparse,logging,subprocess,glob,uuid,shutil,json,pwd
from mysql import mysql_wrapper
logger = logging.getLogger(__name__)
try:
   from rucio.client import Client
   from rucio.common.exception import DataIdentifierNotFound,DataIdentifierAlreadyExists,DuplicateContent
except ImportError,e:
   logger.error('Failed to import Rucio, have you setup Rucio in your environment?')
   raise

DEFAULT_RUCIO_SCOPE='group.phys-gener'
DEFAULT_RUCIO_RESOURCE='ANLASC_DATADISK'
DEFAULT_PATH='/grid/atlas/dq2/ATLASDATADISK/rucio/'
DEFAULT_CONTAINER=None
DEFAULT_USER=os.environ['USER']

MAX_RUCIO_FILES_PER_ADD=500

def main():
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='')
   parser.add_option('-g','--input-glob',dest='input_glob',help='The input file pattern for globbing a file list. Enclose with quotations.')
   parser.add_option('-d','--dataset-name',dest='dataset_name',help='Name of the dataset to upload with.')
   parser.add_option('-c','--container',dest='container_name',help='Name of the container to which to add the dataset. If not set, dataset is not added to a container.',default=DEFAULT_CONTAINER)
   parser.add_option('-s','--rucio-scope',dest='rucio_scope',help='Rucio scope to use during upload [default='+DEFAULT_RUCIO_SCOPE+']',default=DEFAULT_RUCIO_SCOPE)
   parser.add_option('-r','--rucio-rse',dest='rucio_resource',help='Rucio resource to which to upload [default='+DEFAULT_RUCIO_RESOURCE+']',default=DEFAULT_RUCIO_RESOURCE)
   parser.add_option('-p','--path',dest='rucio_path',help='Rucio path on resource to which to upload [default='+DEFAULT_PATH+']',default=DEFAULT_PATH)
   parser.add_option('-u','--use-user',dest='username',help='Execute as username via sudo command. Must execute python using sudo for it to work [default='+DEFAULT_USER+']',default=DEFAULT_USER)
   parser.add_option('--proxyfile',dest='proxyfile',help='If set, changes proxy file from default "/tmp/..." area',default=None)
   parser.add_option('--db-entries',dest='db_entries',help='A comma separated list of DB ids for the Argo job database. If this list is given, the group_identifier will have the text ".uploaded" added to each entry upon the upload successfully completed.',default='')
   parser.add_option('--target-events-total',dest='target_events_total',help='If set, checks if the total number of events have been uploaded to the container in rucio',default=-1,type='int')
   parser.add_option('--events-per-file',dest='events_per_file',help='If target-events-total is set, this is the multiplier. Rucio only returns how many files were uploaded, so is the number multiplied by the number of files returned by rucio to check against target-events-total.',default=5000,type='int')
   options,args = parser.parse_args()

   
   manditory_args = [
                     'input_glob',
                     'dataset_name',
                     'rucio_scope',
                     'rucio_resource',
                     'rucio_path',
                     'username',
                     'proxyfile',
                     'db_entries',
                     'target_events_total',
                     'events_per_file',
                  ]

   for man in manditory_args:
      if options.__dict__[man] is None:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)
   
   file_list = glob.glob(options.input_glob)
   if len(file_list) == 0:
      logger.error(' glob pattern "' + options.input_glob + '" returned no files')
      sys.exit(-1)

   rucio_bulk_upload(
         file_list,
         options.dataset_name,
         options.container_name,
         options.rucio_scope,
         options.rucio_resource,
         options.rucio_path,
         options.username,
         options.proxyfile,
         options.db_entries,
         options.target_events_total,
         options.events_per_file,
      )


def rucio_bulk_upload(
      file_list,
      dataset_name,
      container_name = DEFAULT_CONTAINER,
      rucio_scope = DEFAULT_RUCIO_SCOPE,
      rucio_resource = DEFAULT_RUCIO_RESOURCE,
      rucio_path = DEFAULT_PATH,
      username = DEFAULT_USER,
      proxyfile = None,
      db_entries = '',
      target_events_total = -1,
      events_per_file = 5000,
   ):

   logger.info(' file_list contains ' + str(len(file_list)) + ' files.')
   logger.info('  dataset_name:   ' + dataset_name)
   logger.info('  container_name: ' + container_name)
   logger.info('  rucio_scope:    ' + rucio_scope)
   logger.info('  rucio_resource: ' + rucio_resource)
   logger.info('  rucio_path:     ' + rucio_path)
   logger.info('  username:       ' + username)
   logger.info('  proxyfile:      ' + str(proxyfile))
   logger.info('  db_entries:     ' + db_entries)

   if username != DEFAULT_USER:
      try:
         stage_proxy(DEFAULT_USER,username,proxyfile)
      except:
         raise

   # check for dataset in rucio
   logger.info('check if dataset exists: %s',dataset_name)
   if dataset_exists(dataset_name):
      nfiles = dataset_files(dataset_name)
      logger.warning(' dataset ' + str(dataset_name) + ' already exists in rucio with '+str(nfiles) + ' files compared to input file list which has ' + str(len(file_list)) + '. Exiting.')
      return

   # rucio client
   client = Client()
   
   # the client.add_files_to_dataset can only deal with 
   # 1000 files at a time so we need sub-loops
   logger.info('looping over filelist')
   for n in range(int(len(file_list)/MAX_RUCIO_FILES_PER_ADD)+1):
      logger.info('n = ' + str(n))
      files_done = n*MAX_RUCIO_FILES_PER_ADD
      files_left = len(file_list)-files_done
      if files_left > MAX_RUCIO_FILES_PER_ADD:
         files_left = MAX_RUCIO_FILES_PER_ADD
      logger.info('files_done: %d files_left: %d',files_done,files_left)
      registr_files = []
      for i in range(files_left):
         file_index = i + files_done
         rf = RucioFileAttributes(rucio_scope,file_list[file_index],username)
         rf.copy_file(rucio_scope,rucio_path)
         registr_files.append(rf.get_register_dictionary(rucio_scope))
         logger.info(' i = ' + str(i) + ' file_index = ' + str(file_index))# + ' >> ' + str(rf.get_register_dictionary(rucio_scope)))
      
      
      attachments = {'scope': rucio_scope, 'name': dataset_name, 'rse': rucio_resource, 'dids': registr_files} 
      #logger.info('attachments = %s',json.dumps(attachments))
      #client.add_files_to_dataset(rucio_scope,dataset_name,registr_files,rse=rucio_resource)
      if len(registr_files) > 0:
         try:
            #logger.info('adding %d files to dataset',len(registr_files))
            client.add_files_to_datasets(attachments=[attachments],ignore_duplicate=True) 
         except DataIdentifierNotFound,de:
            logger.info('Adding dataset: ' + rucio_scope + ':' + dataset_name)
            client.add_dataset(scope=rucio_scope,name=dataset_name)
            client.add_files_to_datasets(attachments=[attachments],ignore_duplicate=True)
         except Exception,e:
            logger.exception('Exception raised during add_files_to_datasets. \n attachments = %s' % str(attachments))
            raise
   logger.info('adding dataset to container')
   if container_name is not None:
      try:
         logger.info('Adding container')
         if not client.add_container(rucio_scope,container_name):
            logger.error('add_container command returned False.')
      except DataIdentifierAlreadyExists,e:
         logger.info('container already exists')
      except Exception,e:
         logger.exception('Exception raised during add_container')
         raise

      try:
         logger.info('Attaching dataset to container')
         dids = []
         dids.append({'scope':rucio_scope,'name':dataset_name})
         client.attach_dids(scope=rucio_scope,name=container_name,dids=dids)
      except DuplicateContent,e:
         logger.info('dataset already attached to container')
      except Exception,e:
         logger.exception('Exception raised during attach_dids.')
         raise

   logger.info('upload finished')

   num_files_in_contianer = 0
   try:
      for file_descriptor in client.list_files(rucio_scope,container_name):
         num_files_in_contianer += 1
   except:
      logger.exception('received exception when trying to check file contents of container %s',container_name)

   logger.info('container has %s files',num_files_in_contianer)

   validated = True
   if target_events_total > 0:
      logger.info('validating upload with rucio')
      if (num_files_in_contianer-1) * events_per_file >= target_events_total:
         logger.info('container contents validated. Contains at least events %s events (contents = %s)',target_events_total,(num_files_in_contianer-1)*events_per_file)
      else:
         logger.info('contianer contains less than requested events. Contents = %s, requested = %s',(num_files_in_contianer-1)*events_per_file, target_events_total)
         validated = False


   logger.info('updating group_identifier for ids: '+str(db_entries))
   db_entries = db_entries.split(',')
   logger.info('updating group_identifier for ids: '+str(db_entries))
   if len(db_entries) > 0:
      db = mysql_wrapper.MySQL(db_server_password='1w@ANLHEP')

      for dbid in db_entries:
         logger.info('updating id %s',str(dbid))
         entries,labels = db.select('group_identifier','id=%s'%dbid)
         groupid = entries[0][labels['group_identifier']]
         if '.uploading' in groupid:
            if validated:
               groupid = groupid.replace('.uploading','.uploaded')
            else:
               groupid = groupid.replace('.uploading','.upload_failed')
         else:
            if validated:
               groupid = groupid + '.uploaded'
            else:
               groupid = groupid + '.upload_failed'
         db.update("group_identifier='%s'" % groupid,where_conditions='id=%s' % dbid)

   logger.info(' Done. ')


class RucioFileAttributes:
   def __init__(self,scope,filename,username):
      self.filename  = filename
      self.md5sum    = md5sum(scope,os.path.basename(filename))
      self.uuid      = uuid.uuid1()
      self.bytes     = os.path.getsize(filename)
      self.adler32   = adler32(filename)
      self.username  = username

   def get_file_path(self,scope,path = None):
      filepath = ''
      if path is not None:
         filepath = path

      scope_parts = scope.split('.')
      for part in scope_parts:
         filepath = os.path.join(filepath,part)
      filepath = os.path.join(filepath,
                              str(self.md5sum[0:2]),
                              str(self.md5sum[2:4])
                             )
      return filepath

   def copy_file(self,scope,path = None):

      new_path = self.get_file_path(scope,path)
      new_file = os.path.join(new_path,os.path.basename(self.filename))
      
      if self.username == DEFAULT_USER:
         err = os.system('mkdir -p ' + new_path)
         if err != 0:
            raise Exception('return from "mkdir -p ' + new_path + '" is non-zero: ' + str(err))
         err = os.system('cp ' + self.filename + ' ' + new_file)
         if err != 0:
            raise Exception('return from "cp ' + self.filename + ' ' + new_file + '" is non-zero: ' + str(err))
      else:
         err = os.system('sudo -S -u ' + self.username + ' mkdir -p ' + new_path)
         if err != 0:
            raise Exception('return from "sudo -S -u ' + self.username + ' mkdir -p ' + new_path + '" is non-zero: ' + str(err))
         err = os.system('sudo -S -u ' + self.username + ' cp ' + self.filename + ' ' + new_file)
         if err != 0:
            raise Exception('return from "sudo -S -u ' + self.username + ' cp ' + self.filename + ' ' + new_file + '" is non-zero: ' + str(err))

      if not os.path.exists(new_path):
         raise Exception('Error copying file')

   def get_register_dictionary(self,scope):
      dic = {}
      dic['scope']   = scope
      dic['name']    = os.path.basename(self.filename)
      dic['bytes']   = self.bytes
      dic['adler32'] = self.adler32
      dic['meta']    = {'guid':str(self.uuid)}
      return dic

   def __str__(self):
      fparts = self.filename.split('/')
      tmp  = ' filename: ' + fparts[0]
      if len(fparts) > 1:
         for i in range(1,len(fparts)):
            tmp += '/\n             ' + fparts[i]
      tmp += '\n'
      tmp += ' md5sum:   ' + self.md5sum + '\n'
      tmp += ' uuid:     ' + str(self.uuid) + '\n'
      tmp += ' bytes:    ' + str(self.bytes) + '\n'
      tmp += ' adler32:   0x' + self.adler32 + '\n'
      return tmp


def md5sum(scope,filename):
   cmd = 'echo -n "' + scope + ':' + filename + '" | md5sum'
   p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
   p.wait()
   if p.returncode != 0:
      raise Exception('md5sum returned non-zero code: '+ str(p.returncode))
   stdout,stderr = p.communicate()

   return stdout.split()[0]

def adler32(filename):
   exe = 'xrdadler32'
   cmd = exe + ' ' + filename
   p = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
   stdout,stderr = p.communicate()
   if p.returncode != 0:
      raise Exception('adler32 returned non-zero code: %s \n stdout: %s \n stderr: %s \n cmd: %s' % (str(p.returncode),stdout,stderr,cmd))
   

   return stdout.split()[0]


class ProxyDoesNotExist(Exception): pass
class ErrorCopyingProxy(Exception): pass
class ErrorChownProxy(Exception): pass
def stage_proxy(from_username,to_username,proxyfile = None):

   
   if proxy_exists(proxyfile):
      if proxyfile is None:
         from_proxy = 'x509up_u' + str(pwd.getpwnam(from_username).pw_uid)
      else:
         from_proxy = proxyfile
      to_proxy   = 'x509up_u' + str(pwd.getpwnam(to_username).pw_uid)

      cmd = 'sudo cp ' + os.path.join('/tmp',from_proxy) + ' ' + os.path.join('/tmp',to_proxy)
      logger.info('staging proxy: ' + cmd)
      p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
      p.wait()
      stdout,stderr = p.communicate()
      if p.returncode != 0:
         raise ErrorCopyingProxy('error copying proxy of user ' + from_username + ' to user ' 
            + to_username + ' stdout: ' + stdout)

      cmd = 'sudo chown ' + to_username + ' ' + os.path.join('/tmp',to_proxy)
      p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
      p.wait()
      stdout,stderr = p.communicate()
      if p.returncode != 0:
         raise ErrorChownProxy('error changing ownership of proxy ' + os.path.join('/tmp',to_proxy) 
            + ' to user ' + to_username + ' stdout: ' + stdout)

   else:
      raise ProxyDoesNotExist('proxy does not exist')

def proxy_exists(proxyfile = None):
   cmd = 'grid-proxy-info -exists'
   if proxyfile is not None:
      cmd = 'X509_USER_PROXY=%s %s' % (proxyfile,cmd)
   p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
   p.wait()
   if p.returncode == 0:
      return True
   return False

def dataset_exists(dataset_name):
   return False
   cmd = 'rucio list-dids ' + dataset_name
   logger.info('checking dataset exists: %s',cmd)
   p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
   p.wait()
   if p.returncode == 0:
      return True
   return False

def dataset_files(dataset_name):
   cmd = 'rucio list-files ' + dataset_name
   p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
   stdout,stderr = p.communicate()
   if p.returncode == 0:
      start_index = stdout.rfind('Total files :')
      end_index = stdout.rfind('Total size :')
      if start_index > 0 and end_index > 0:
         return int(stdout[start_index+13:end_index-1])
   return 0

if __name__ == "__main__":
   main()
