#!/usr/bin/env python

import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

import optparse,os,sys,time,shutil

from ArgoJobs.PythiaArgoJob import PythiaArgoJob

sys.path.append('/users/hpcusers/balsam/balsam_deploy/common_core')
from MessageInterface import MessageInterface
from ArgoJob import ArgoJob

def main():
   # parse command line
   parser = optparse.OptionParser(description='job payload for running Pythia8')
   parser.add_option('-i','--input-file',dest='input_file',help='filename of the input configuration card for pythia')
   parser.add_option('-n','--num-evts',dest='numevts',help='number of events to generate',type='int')
   parser.add_option('-s','--site',dest='site',help='Balsam site name on which to run the event generation')
   parser.add_option('-o','--num-nodes',dest='numnodes',help='number of nodes to use on destination machine',type='int')
   parser.add_option('-c','--cpus-per-node',dest='cpus_per_node',help='number of CPUs per node to use on destination machine',type='int')
   parser.add_option('-t','--wall-minutes',dest='wall_minutes',help='number of wall minutes your job will be assigned in the scheduler',type='int')
   parser.add_option('-x','--no-submit',dest='submit',help='do not submit the message to ARGO. For testing purposes.',action='store_false',default=True)
   parser.add_option('-d','--dev',dest='dev',help='use development servers/rabbit',action='store_true',default=False)
   options,args = parser.parse_args()
   
   if options.numevts is None:
      parser.error('Must define the number of events')
   elif options.site is None:
      parser.error('Must define the destination site')
   elif options.input_file is None:
      parser.error('Must define the input configuration filename')
   elif not options.numnodes:
      parser.error('Must define the number of nodes to use')
   elif not options.cpus_per_node:
      parser.error('Must define the number of CPUs per node to use')
   elif not options.wall_minutes:
      parser.error('Must define the wall minutes')
   
   submit_pythia(
                  options.input_file,
                  options.numevts,
                  options.site,
                  options.numnodes,
                  options.cpus_per_node,
                  options.wall_minutes,
                  options.submit,
                  options.dev
                )



def submit_pythia(
                  input_filename,
                  nevts,
                  site,
                  nodes,
                  ranks_per_node,
                  wall_minutes,
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
   
   job = PythiaArgoJob()
   job.input_filename   = input_filename
   job.working_path     = RUNPATH
   job.nevts            = nevts
   job.nodes            = nodes
   job.ranks_per_node   = ranks_per_node
   job.site             = site
   job.input_url        = input_url
   job.output_url       = output_url
   job.wall_minutes     = wall_minutes

   argojob = job.get_argo_job()
   argojob.email_address = 'jchilders@anl.gov,turam@anl.gov'
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



