#!/usr/bin/env python
import os,sys,optparse,logging,subprocess,glob,uuid,shutil,json,pwd
logger = logging.getLogger(__name__)


DEFAULT_RUCIO_SCOPE='group.phys-gener'
DEFAULT_RUCIO_RESOURCE='ANLASC_SCRATCHDISK'
DEFAULT_PATH='/grid/atlas/dq2/ATLASSCRATCHDISK/rucio/'
DEFAULT_CONTAINER=None
DEFAULT_USER=os.environ['USER']

MAX_RUCIO_FILES_PER_ADD = 499

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
   options,args = parser.parse_args()

   
   manditory_args = [
                     'input_glob',
                     'dataset_name',
                     'rucio_scope',
                     'rucio_resource',
                     'rucio_path',
                     'username',
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
      )


def rucio_bulk_upload(
      file_list,
      dataset_name,
      container_name = DEFAULT_CONTAINER,
      rucio_scope = DEFAULT_RUCIO_SCOPE,
      rucio_resource = DEFAULT_RUCIO_RESOURCE,
      rucio_path = DEFAULT_PATH,
      username = DEFAULT_USER,
   ):

   logger.info(' file_list contains ' + str(len(file_list)) + ' files.')
   logger.info('  dataset_name:   ' + dataset_name)
   logger.info('  container_name: ' + container_name)
   logger.info('  rucio_scope:    ' + rucio_scope)
   logger.info('  rucio_resource: ' + rucio_resource)
   logger.info('  rucio_path:     ' + rucio_path)
   logger.info('  username:       ' + username)

   if username != DEFAULT_USER:
      try:
         stage_proxy(DEFAULT_USER,username)
      except:
         raise

   # check for dataset in rucio
   # if dataset_exists(dataset_name):
   #    nfiles = dataset_files(dataset_name)
   #    logger.warning(' dataset ' + str(dataset_name) + ' already exists in rucio with '+str(nfiles) + ' files compared to input file list which has ' + str(len(file_list)) + '. Exiting.')
   #    return

   # rucio client
   # client = Client()
   
   # the client.add_files_to_dataset can only deal with 
   # 1000 files at a time so we need sub-loops
   for n in range(int(len(file_list)/MAX_RUCIO_FILES_PER_ADD)+1):
      logger.info('n = ' + str(n))
      files_done = n*MAX_RUCIO_FILES_PER_ADD
      files_left = len(file_list)-files_done
      if files_left > MAX_RUCIO_FILES_PER_ADD:
         files_left = MAX_RUCIO_FILES_PER_ADD
      #registr_files = []
      for i in range(files_left):
         file_index = i + files_done
         rf = RucioFileAttributes(rucio_scope,file_list[file_index],username)
         logger.info(' i = ' + str(i) + ' file_index = ' + str(file_index))
         rf.copy_file(rucio_scope,rucio_path)
         #registr_files.append(rf.get_register_dictionary(rucio_scope))

      
      #attachments = {'scope': rucio_scope, 'name': dataset_name, 'rse': rucio_resource, 'dids': registr_files} 
      #print json.dumps(attachments)
      #client.add_files_to_dataset(rucio_scope,dataset_name,registr_files,rse=rucio_resource)
      # try:
      #    client.add_files_to_datasets(attachments=[attachments],ignore_duplicate=True) 
      # except DataIdentifierNotFound,de:
      #    logger.info('Adding dataset: ' + rucio_scope + ':' + dataset_name)
      #    client.add_dataset(scope=rucio_scope,name=dataset_name)
      #    client.add_files_to_datasets(attachments=[attachments],ignore_duplicate=True)
      # except Exception,e:
      #    logger.exception('Exception raised during add_files_to_datasets.')
      #    raise

   # if container_name is not None:
   #    try:
   #       logger.info('Adding container')
   #       if not client.add_container(rucio_scope,container_name):
   #          logger.error('add_container command returned False.')
   #    except DataIdentifierAlreadyExists,e:
   #       logger.info('container already exists')
   #    except Exception,e:
   #       logger.exception('Exception raised during add_container')
   #       raise

   #    try:
   #       logger.info('Attaching dataset to container')
   #       dids = []
   #       dids.append({'scope':rucio_scope,'name':dataset_name})
   #       client.attach_dids(scope=rucio_scope,name=container_name,dids=dids)
   #    except DuplicateContent,e:
   #       logger.info('dataset already attached to container')
   #    except Exception,e:
   #       logger.exception('Exception raised during attach_dids.')
   #       raise

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

      logger.info('trying to copy: \n' + self.filename + '\n to ' + new_file )

      if os.path.exists(self.filename) and not os.path.exists(new_file):
         
         if self.username == DEFAULT_USER:
            os.system('mkdir -p ' + new_path)
            os.system('cp ' + self.filename + ' ' + new_file)
         else:
            os.system('sudo -S -u ' + self.username + ' mkdir -p ' + new_path)
            os.system('sudo -S -u ' + self.username + ' cp ' + self.filename + ' ' + new_file)
      elif not os.path.exists(self.filename):
         logger.error(' source file doesn not exist.')
      elif os.path.exists(new_file):
         logger.error(' destination already exists.')

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
   cmd = 'xrdadler32 ' + filename
   p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
   p.wait()
   if p.returncode != 0:
      raise Exception('adler32 returned non-zero code: '+ str(p.returncode))
   stdout,stderr = p.communicate()

   return stdout.split()[0]


class ProxyDoesNotExist(Exception): pass
class ErrorCopyingProxy(Exception): pass
class ErrorChownProxy(Exception): pass
def stage_proxy(from_username,to_username):

   
   if proxy_exists():
      from_proxy = 'x509up_u' + str(pwd.getpwnam(from_username).pw_uid)
      to_proxy   = 'x509up_u' + str(pwd.getpwnam(to_username).pw_uid)

      cmd = 'sudo cp ' + os.path.join('/tmp',from_proxy) + ' ' + os.path.join('/tmp',to_proxy)
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

def proxy_exists():
   cmd = 'grid-proxy-info -exists'
   p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
   p.wait()
   if p.returncode == 0:
      return True
   return False

def dataset_exists(dataset_name):
   cmd = 'rucio list-dids ' + dataset_name
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
