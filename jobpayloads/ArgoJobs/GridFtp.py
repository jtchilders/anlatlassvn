import sys,os,subprocess,shlex,logging
#from django.conf import settings

class settings:
   path = '/usr/bin/'
   GRIDFTP_PROXY_INFO = path + 'grid-proxy-info'
   GRIDFTP_PROXY_INIT = path + 'grid-proxy-init'
   GRIDFTP_GLOBUS_URL_COPY = path + 'globus-url-copy'

logger = logging.getLogger(__name__)

# job has completed, copy files from gridftp server
#logger.debug('setup X509 certificates for gridftp calls')
#os.environ['X509_USER_CERT'] = '/users/hpcusers/balsam/gridsecurity/jchilders/xrootdsrv-cert.pem'
#os.environ['X509_USER_KEY']  = '/users/hpcusers/balsam/gridsecurity/jchilders/xrootdsrv-key.pem'

def globus_url_copy(from_path,to_path):

   p = subprocess.Popen([settings.GRIDFTP_PROXY_INFO,'-exists'])
   p.wait()
   if p.returncode is not 0: # valid proxy does not exist so create one
      p = subprocess.Popen([settings.GRIDFTP_PROXY_INIT,'-verify','-debug','-bits','2048','-valid','96:00'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      stdout = ''
      stderr = ''
      for line in p.stdout:
         logger.info(line[0:-1])
         stdout += line
      for line in p.stderr:
         logger.info(line[0:-1])
         stderr += line
      p.wait()
      if p.returncode is not 0:
         logger.warning('grid-proxy-init failed: stdout = \n' + stdout + '\n stderr: ' + stderr + '\n')
      else:
         logger.debug(' grid-proxy-init stdout: \n' + stdout + '\n stderr: \n' + stderr + '\n')


   from_path = str(from_path)
   to_path = str(to_path)
   cmd = settings.GRIDFTP_GLOBUS_URL_COPY + ' -cd  -r -nodcau ' + from_path + ' ' + to_path
   #logger.debug(' GridFtp, from_path = ' + from_path )
   #logger.debug(' GridFtp, to_path   = ' + to_path )

   logger.info(' transferring files: ' + cmd )
   try:
      proc = subprocess.Popen(shlex.split(cmd),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   except:
      logger.exception(' Error while performing GridFTP transfer: ' + str(sys.exc_info()[1]) )
      raise
   stdout = ''
   stderr = ''
   for line in proc.stdout:
      logger.info(line[0:-1])
      stdout += line
   for line in proc.stderr:
      logger.info(line[0:-1])
      stderr += line
   p.wait()
   return stdout,stderr


