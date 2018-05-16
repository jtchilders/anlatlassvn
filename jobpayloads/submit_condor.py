#!/usr/bin/env python

import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

import optparse,os,sys,time,shutil

from ArgoJobs.CondorArgoJob import CondorArgoJob

sys.path.append('/users/hpcusers/balsam/balsam_deploy/common_core')
from MessageInterface import MessageInterface
from ArgoJob import ArgoJob

def main():
   # parse command line
   parser = optparse.OptionParser(description='job payload for running simple condor and condor dagman jobs')
   parser.add_option('-a','--condor-job',dest='condor_job_file',help='filename of the input condor job')
   parser.add_option('-b','--condor_dagman',dest='condor_dagman_file',help='filename of the input condor dagman job')
   parser.add_option('-i','--input-files',dest='input_files',help='comma separated list of files that are needed for the job')
   parser.add_option('-j','--output-files',dest='output_files',help='comma separated list of files that will be output from the job')
   #parser.add_option('-o','--num-nodes',dest='numnodes',help='number of nodes to use on destination machine',type='int')
   #parser.add_option('-c','--cpus-per-node',dest='cpus_per_node',help='number of CPUs per node to use on destination machine',type='int')
   #parser.add_option('-t','--wall-minutes',dest='wall_minutes',help='number of wall minutes your job will be assigned in the scheduler',type='int')
   parser.add_option('-x','--no-submit',dest='submit',help='do not submit the message to ARGO. For testing purposes.',action='store_false',default=True)
   parser.add_option('-d','--dev',dest='dev',help='use development servers/rabbit',action='store_true',default=False)
   options,args = parser.parse_args()
   
   if options.condor_job_file is None and options.condor_dagman_file is None:
      parser.error('Must define an input job file or DAGMAN job file')
   
   input_files = None
   if options.input_files:
      input_files = options.input_files.split(',')
   else:
      input_files = []
   output_files = None
   if options.output_files:
      output_files = options.output_files.split(',')
   else:
      output_files = []

   submit_condor(
                  options.condor_job_file,
                  options.condor_dagman_file,
                  input_files,
                  output_files,
                  options.submit,
                  options.dev
                )



def submit_condor(
                  condor_job_file,
                  condor_dagman_file,
                  input_files = [],
                  output_files = [],
                  submit = True,
                  dev = False,
                 ):

   # Assigns a unique task number.  This is a sixteen digit number.
   
   # For non-PanDA submission the top ten digits represent
   # the time of submission, with the bottom six being random
   # For PanDA submission, the top ten digits are the PanDA
   # job ID and the next six are zeroes.
   pandaID = os.environ.get('PandaID', None)
  
   if ( pandaID != None):
      taskID = pandaID + '000000'
   else:
      taskID=str(int(time.time()*1000000));

   # Find out who is submitting this job and place it in 
   #   the variable 'user'. For some reason prodUserID
   #   is not set.

   user = os.environ.get('USER','nobody')
   if(user == 'apf'):
      user= os.environ.get('prodUserID','nobody')
   jobID = taskID + '0'
   
   # create gridFTP URLs 
   grid_ftp_server = 'atlasgridftp02.hep.anl.gov'
   grid_base_path  = '/grid/atlas/hpc/argo/jobs/' + str(jobID)
   input_url = 'gsiftp://' + grid_ftp_server + grid_base_path
   output_url = input_url
   
   TOP_PATH = os.getcwd() # directory in which script was run
   RUNPATH = os.path.join(TOP_PATH,str(jobID)) # directory in which to store files
   if not os.path.exists(RUNPATH):
      logger.info('creating run path: ' + RUNPATH)
      os.makedirs(RUNPATH) # make directories recursively like 'mkdir -p'
   
   condor_job = CondorArgoJob()
   condor_job.condor_job_file     = condor_job_file
   condor_job.condor_dagman_file  = condor_dagman_file
   condor_job.working_path        = RUNPATH
   condor_job.input_files         = input_files
   condor_job.output_files        = output_files
   condor_job.input_url           = input_url
   condor_job.output_url          = output_url

   if dev:
      condor_job.site = 'argo_cluster_dev'

   argojob = condor_job.get_argo_job()
   argojob.email_address = 'jchilders@anl.gov'
   job_txt = argojob.serialize()

   print job_txt

   if submit:
      mi = MessageInterface()
      mi.host = 'atlasgridftp02.hep.anl.gov'
      mi.port = 5671
      mi.ssl_cert = os.environ['X509_USER_CERT'] #'/users/hpcusers/balsam/gridsecurity/jchilders/xrootdsrv-cert.pem'
      mi.ssl_key  = os.environ['X509_USER_KEY'] #'/users/hpcusers/balsam/gridsecurity/jchilders/xrootdsrv-key.pem'
      mi.ssl_ca_certs = os.environ['X509_CA_CERTS'] #'/users/hpcusers/balsam/gridsecurity/jchilders/cacerts.pem'
      mi.exchange_name = 'argo_users'
      if dev:
         mi.exchange_name = 'argo_users_dev'
      logger.info('opening connection')
      mi.open_blocking_connection()
      routing_key = 'argo_job'
      if dev:
         routing_key = 'argo_job_dev'
      logger.info(' sending msg ')
      mi.send_msg(argojob.serialize(),routing_key)
      logger.info( ' done sending')
      mi.close()
      logger.info(' closing connection')
   else:
      logger.info(' not submitting job ')


if __name__ == '__main__':
   main()



