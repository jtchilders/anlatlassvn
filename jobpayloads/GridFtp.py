import sys,os,subprocess,shlex,logging
logger = logging.getLogger(__name__)

# job has completed, copy files from gridftp server
logger.debug('setup X509 certificates for gridftp calls')
os.environ['X509_USER_CERT'] = '/users/hpcusers/balsam/gridsecurity/jchilders/xrootdsrv-cert.pem'
os.environ['X509_USER_KEY']  = '/users/hpcusers/balsam/gridsecurity/jchilders/xrootdsrv-key.pem'

p = subprocess.Popen(['grid-proxy-info','-exists'])
p.wait()
if p.returncode is not 0: # valid proxy does not exist so create one
   p = subprocess.Popen(['grid-proxy-init','-debug'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   stdout,stderr = p.communicate()
   if p.returncode is not 0:
      logger.warning('grid-proxy-init failed: stdout = \n' + stdout + '\n stderr: ' + stderr + '\n')
   else:
      logger.debug(' grid-proxy-init stdout: \n' + stdout + '\n stderr: \n' + stderr + '\n')
   



def globus_url_copy(from_path,to_path):
   cmd = 'globus-url-copy -cd  -r -nodcau ' + from_path + ' ' + to_path
   logger.debug(' transfering files: ' + cmd)
   proc = subprocess.Popen(shlex.split(cmd),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   stdout,stderr = proc.communicate()
   logger.debug(' transfer stdout: \n' + stdout + '\n stderr: \n' + stderr + '\n')
   return stdout,stderr


