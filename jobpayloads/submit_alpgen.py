#!/usr/bin/env python
import os,sys,logging,optparse,time,logging,json,shutil,glob,random,datetime

sys.path.append('/users/hpcusers/svn/generators/alpgen/v214/usercode/python')
from ArgoJobs.AlpgenArgoJob import AlpgenArgoJob
from ArgoJob import ArgoJob

sys.path.append('/users/hpcusers/svn/generators/alpgen/v214/usercode/python')
import AlpgenInputFile

logger = logging.getLogger(__name__)

JOB_DEPOT = 'jobs_alpgen'

random.seed(datetime.datetime.now())
def get_rand():
   return random.randrange(2**16) # alpgen requires iseeds to be 5-digit integer

def main():
   
   logging.basicConfig(level=logging.INFO)
   parser = optparse.OptionParser(description='submit alpgen job to ARGO')
   parser.add_option('-e','--evts-per-iter',dest='evts_per_iter',help='number of events per warmup iteration',type='int')
   parser.add_option('-i','--num-iter',dest='numiters',help='number of iterations for the warmup',type='int')
   parser.add_option('-w','--warmup-weighted',dest='num_warmup',help='number of event to in the warmup, after the iterations complete',type='int')
   parser.add_option('-n','--num-weighted',dest='num_weighted',help='number of weighted events to generate.',type='int')
   parser.add_option('-p','--process',dest='process',help='define the process to generate, 2Q,4Q,hjet,top,wjet,zjet,Njet,etc.')
   parser.add_option('-o','--num-nodes',dest='numnodes',help='number of nodes to use on destination machine',type='int')
   parser.add_option('-c','--cpus-per-node',dest='cpus_per_node',help='number of CPUs per node to use on destination machine',type='int')
   parser.add_option('-a','--alpgen-input',dest='alpgen_input_file',help='The AlpGen input file which carries all the options for this generation job')
   parser.add_option('-f','--pdf-filename',dest='pdf_filename',help='The PDF Filename for Alpgen, default="./cteq6l1.tbl"',default='cteq6l1.tbl')
   parser.add_option('-t','--wall-time',dest='walltime',help='The wall time to submit to the queue in minutes.',type='int')
   parser.add_option('-s','--site',dest='site',help='Balsam site name on which to run the event generation')
   parser.add_option('-d','--dev',dest='dev',help='Run in development mode, means warmup is sent to argo_cluster_dev',action='store_true',default=False)
   parser.add_option('-r','--regen',dest='regen_path',help='Regenerate the events from a previous warmup, the argument is the path from which to copy the grid1/2 files. The input file should still be provided using the -a option.')
   parser.add_option('-x','--no-submit',dest='submit',help='do not submit the message to ARGO. For testing purposes.',action='store_false',default=True)
   parser.add_option('','--iseed1',dest='iseed1',help='override random number, iseed1, in alpgen input file',type='int')
   parser.add_option('','--iseed2',dest='iseed2',help='override random number, iseed2, in alpgen input file',type='int')
   parser.add_option('','--iseed3',dest='iseed3',help='override random number, iseed3, in alpgen input file',type='int')
   parser.add_option('','--iseed4',dest='iseed4',help='override random number, iseed4, in alpgen input file',type='int')
   parser.add_option('-g','--group-id',dest='group_identifier',help='This identifier is added to the ARGO database and can be used to label jobs for reference later.')
   parser.add_option('-q','--status-queue',dest='enable_status_queue',help='Enable the setting of the message queue parameter in the ArgoJob, which means ARGO will not send message updates for this job to the queue with its job ID.',action='store_true',default=False)
   parser.add_option('--njobs',dest='njobs',help='Use this option to submit more than one job of this type. The random number will be updated automatically for each job.',type='int',default=1)
   parser.add_option('--norandomize',dest='randomize',help='Turns off creation of random seeds for each job',action='store_false',default=True)

   options,args = parser.parse_args()


   if options.numiters is None and options.regen_path is None:
      parser.error('Must define the number of warmup iterations')
   if options.process is None:
      parser.error('Must define the process to generate')
   if options.numnodes is None:
      parser.error('Must define the number of nodes to use')
   if options.cpus_per_node is None:
      parser.error('Must define the number of CPUs per node to use')
   if options.evts_per_iter is None and options.regen_path is None:
      parser.error('Must define the number of events per warmup iteration')
   if options.num_weighted is None:
      parser.error('Must define the number of weighted events to produce')
   if options.num_warmup is None and options.regen_path is None:
      parser.error('Must define the number of weighted events to produce in the warmup step.')
   if options.alpgen_input_file is None and options.regen_path is None:
      parser.error('Must define the AlpGen input file')
   if options.walltime is None:
      parser.error('Must specify a wall time')
   if options.pdf_filename is None:
      parser.error('Must specify a PDF Filename')
   if options.walltime is None:
      parser.error('Must specify a wall time')
   if options.site is None:
      parser.error('Must specify a site')
   
   
   if options.site == 'mira':
      if( options.numnodes != 2**9
          and options.numnodes != 2**10
          and options.numnodes != 2**11
          and options.numnodes != 2**12
          and options.numnodes != 2**13
          and options.numnodes != 2**14
          and options.numnodes != 2**15
          and options.numnodes != (2**15 + 2**14)
        ):
        parser.error(' You selected Mira as a site and are using a number of nodes, ' + str(options.numnodes) + ', that is not a power of 2 (512,1024,2048,4096,8192...)')

   regen_path = options.regen_path
   for i in range(options.njobs):
      regen_path = submit_alpgen(
                 options.evts_per_iter,
                 options.numiters,
                 options.num_warmup,
                 options.num_weighted,
                 options.process,
                 options.numnodes,
                 options.cpus_per_node,
                 options.alpgen_input_file,
                 options.pdf_filename,
                 options.walltime,
                 options.site,
                 options.dev,
                 regen_path,
                 options.submit,
                 options.iseed1,
                 options.iseed2,
                 options.iseed3,
                 options.iseed4,
                 options.group_identifier,
                 options.enable_status_queue,
                 options.randomize,
                )


def submit_alpgen(
                 evts_per_iter,
                 numiters,
                 num_warmup,
                 num_weighted,
                 process,
                 numnodes,
                 cpus_per_node,
                 alpgen_input_file,
                 pdf_filename,
                 walltime,
                 site,
                 dev                   = False,
                 regen_path            = None,
                 submit                = True,
                 iseed1                = None,
                 iseed2                = None,
                 iseed3                = None,
                 iseed4                = None,
                 group_identifier      = None,
                 enable_status_queue   = False,
                 randomize             = False,
                 ):

   taskID=str(int(time.time()*1000000))

   if randomize and (iseed1 is not None or iseed2 is not None or iseed3 is not None or iseed4 is not None):
      raise Exception('randomize option set AND iseed values specified... these are mutually exclusive options')
   
   user = os.environ.get('USER','nobody')
   
   jobID = create_job_path() 

   logger.info('JobID: ' + str(jobID))

   TOP_PATH = os.getcwd() # directory in which script was run
   RUNPATH = os.path.join(TOP_PATH,JOB_DEPOT,str(jobID)) # directory in which to store files
   logger.info('Run path: ' + RUNPATH)
   if not os.path.exists(RUNPATH):
      logger.info('creating run path: ' + RUNPATH)
      os.makedirs(RUNPATH) # make directories recursively like 'mkdir -p'
   
   if regen_path is not None:
      logger.info('regenerating from this path: ' + str(regen_path))
      grid1 = glob.glob(regen_path + '/*.grid1.presubmit')[0]
      grid2 = glob.glob(regen_path + '/*.grid2.presubmit')[0]
      logger.info('    copying files: \n' + grid1 + '\n' + grid2)
      shutil.copy(grid1,RUNPATH+'/'+os.path.basename(grid1[0:grid1.find('.presubmit')]))
      shutil.copy(grid2,RUNPATH+'/'+os.path.basename(grid2[0:grid2.find('.presubmit')]))
      # also make copies of the grid files for reference
      shutil.copy(grid1,RUNPATH+'/alpout.grid1')
      shutil.copy(grid1,RUNPATH+'/alpout.grid1.presubmit')
      shutil.copy(grid2,RUNPATH+'/alpout.grid2')
      shutil.copy(grid2,RUNPATH+'/alpout.grid2.presubmit')
      
      # update old input file with incremented iseed1
      alpgenInputFile = AlpgenInputFile.AlpgenInputFile()
      alpgenInputFile.read(os.path.join(regen_path,'alpout.input.1'))
      if 'iseed1' in alpgenInputFile.options:
         value = int(alpgenInputFile.options['iseed1'].value)
         if iseed1 is not None:
            new_value = iseed1
         elif randomize:
            new_value = get_rand()
            iseed1 = new_value # set this so it is picked up later
         else:
            new_value = value + 1

         logger.info('   updating iseed1 from ' + str(value) + ' to ' + str(new_value))
         alpgenInputFile.options['iseed1'].value = str(new_value)
      tmp_input_filename = '/tmp/alpgen.input.pid'+str(os.getpid())
      alpgenInputFile.write(tmp_input_filename)
      alpgen_input_file = tmp_input_filename

      # grab old warmup values for consistency
      alpgenInputFile = AlpgenInputFile.AlpgenInputFile()
      alpgenInputFile.read(os.path.join(regen_path,'alpout.input.0'))
      evts_per_iter = int(alpgenInputFile.nevt)
      numiters = int(alpgenInputFile.nitr)
      num_warmup = int(alpgenInputFile.last_nevt)
      alpgenInputFile = None


      # copy old condor logs for prosperity
      condor_files = glob.glob(os.path.join(regen_path,'condor_*.txt*'))
      for condor_file in condor_files:
         shutil.copy(condor_file,RUNPATH+'/')




   # create input/output gridftp urls
   grid_ftp_server = 'atlasgridftp02.hep.anl.gov'
   grid_base_path  = '/grid/atlas/hpc/argo/jobs/' + str(jobID)
   input_url = 'gsiftp://' + grid_ftp_server + grid_base_path
   output_url = input_url

   # don't want "None" listed for unspecified options, want 0
   if evts_per_iter is None: evts_per_iter = 0
   if numiters is None: numiters = 0
   if num_warmup is None: num_warmup = 0

   job = AlpgenArgoJob()
   job.process                         = process
   job.input_filename                  = alpgen_input_file
   job.warmup_phase0_number_events     = evts_per_iter
   job.warmup_phase0_number_iterations = numiters
   job.warmup_phase1_number_events     = num_warmup
   job.warmup_wall_minutes             = 60
   job.evtgen_phase0_number_events     = 0
   job.evtgen_phase0_number_iterations = 0
   job.evtgen_phase1_number_events     = num_weighted
   job.evtgen_nodes                    = numnodes
   job.evtgen_processes_per_node       = cpus_per_node
   job.evtgen_wall_minutes             = walltime
   job.working_path                    = RUNPATH
   job.input_url                       = input_url
   job.output_url                      = output_url
   job.pdf_filename                    = pdf_filename
   job.username                        = os.environ['USER']
   job.group_identifier                = group_identifier
   if dev:
      job.warmup_site                  = 'argo_cluster_dev'
   job.evtgen_site                     = site
   if iseed1 is not None:
      job.input_iseed1 = iseed1
   elif randomize:
      job.input_iseed1 = get_rand()
   if iseed2 is not None:
      job.input_iseed2 = iseed2
   elif randomize:
      job.input_iseed2 = get_rand()
   if iseed3 is not None:
      job.input_iseed3 = iseed3
   elif randomize:
      job.input_iseed3 = get_rand()
   if iseed4 is not None:
      job.input_iseed4 = iseed4
   elif randomize:
      job.input_iseed4 = get_rand()
   
   argojob = job.get_argo_job()
   if argojob is None:
      logger.error(' error getting argo job ')
      return
   argojob.email_address = 'jchilders@anl.gov,turam@anl.gov'
   # if this is a regeneration using a previously made grid1/2
   # remove the balsam job for the warmup and only do the event gen
   if regen_path is not None:
      argojob.jobs.pop(0)

   # load message interface, but choose the development version if we are developing:
   if dev:
      sys.path.append('/users/hpcusers/balsam_dev/balsam_deploy/common_core')
   else:
      sys.path.append('/users/hpcusers/balsam_production/balsam_deploy/common_core')
   from MessageInterface import MessageInterface

   if enable_status_queue:
      argojob.job_status_routing_key = 'test_job_status' #'status_' + jobID
      logger.info('setting up job status queue: ' + argojob.job_status_routing_key)
      mi = MessageInterface()
      mi.host = 'atlasgridftp02.hep.anl.gov'
      mi.port = 5671
      mi.ssl_cert = os.environ['X509_USER_CERT'] #'/users/hpcusers/balsam/gridsecurity/jchilders/xrootdsrv-cert.pem'
      mi.ssl_key  = os.environ['X509_USER_KEY'] #'/users/hpcusers/balsam/gridsecurity/jchilders/xrootdsrv-key.pem'
      mi.ssl_ca_certs = os.environ['X509_CA_CERTS'] #'/users/hpcusers/balsam/gridsecurity/jchilders/cacerts.pem'
      mi.exchange_name = 'argo_users'
      if dev:
         mi.exchange_name = 'argo_users_dev'
      mi.open_blocking_connection()
      mi.create_queue(argojob.job_status_routing_key,argojob.job_status_routing_key)
      mi.close()

   
   job_txt = argojob.serialize()

   #print job_txt

   x = '''
   print ' '

   new_job = ArgoJob()
   new_job.deserialize(job_txt)

   print new_job.serialize()
   
   if job_txt.find(new_job.serialize()) >= 0:
      print 'identical'

   return '''
   
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
      print 'opening connection'
      mi.open_blocking_connection()
      routing_key = 'argo_job'
      if dev:
         routing_key = 'argo_job_dev'
      print ' sending msg '
      mi.send_msg(argojob.serialize(),routing_key)
      print ' done sending'
      mi.close()
      print ' closing connection'
   else:
      logger.info(' not submitting job ')
      logger.info(job_txt)

   return RUNPATH

def create_job_path():
   job_id = int(time.time()*1e6)
   if not os.path.exists(JOB_DEPOT):
      os.mkdir(JOB_DEPOT)
   # check for existing folders
   while os.path.exists(os.path.join(JOB_DEPOT,str(job_id))):
      job_id = int(time.time()*1e6)
   os.mkdir(os.path.join(JOB_DEPOT,str(job_id)))
   return job_id


if __name__ == '__main__':
   main()
