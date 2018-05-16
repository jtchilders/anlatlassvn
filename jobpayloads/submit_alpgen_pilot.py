#!/usr/bin/env python
import os,sys,logging,optparse,time,logging,json,shutil,glob

sys.path.append('/users/hpcusers/svn/generators/alpgen/v214/usercode/python')
from ArgoJobs.AlpgenArgoJob import PilotAlpgenArgoJob
from ArgoJob import ArgoJob

sys.path.append('/users/hpcusers/svn/generators/alpgen/v214/usercode/python')
import AlpgenInputFile

logger = logging.getLogger(__name__)

JOB_DEPOT = 'jobs_alpgen'

def main():
   
   logging.basicConfig(level=logging.INFO)
   parser = optparse.OptionParser(description='submit alpgen pilot job to ARGO')
   parser.add_option('-c','--ecmEnergy',dest='ecmEnergy',help='Center of mass energy.',type='int')
   parser.add_option('-r','--runNumber',dest='runNumber',help='Job run number.')
   parser.add_option('-j','--jobConfig',dest='jobConfig',help='A comma-separated list of job configuration script files')
   parser.add_option('-e','--evgenJobOpts',dest='evgenJobOpts',help='Download and install the EvgenJobOpts tarball with the given name.')
   parser.add_option('-o','--outputEVNTFile',dest='outputEVNTFile',help='POOL file into which generated events will be written.',type='int')
   parser.add_option('-a','--AMITag',dest='AMITag',help='AMI tag from which this job was defined - this option simply writes the relevant AMI tag value into the output metadata, it does not configure the job.')
   parser.add_option('-s','--steering',dest='steering',help="Steer the transform by manipulating the execution graph before the execution path is calculated. Format is substep:{in,out}{+-}DATA,{in,out}{+-}DATA,... to modify the substep's input/output by adding/removing a data type. e.g. RAWtoESD:in-RDO,in+RDO_TRIG would remove RDO and add RDO_TRIG to the list of valid input datatypes for the RAWtoESD substep..")
   parser.add_option('-g','--gridgen-site',dest='gridgen_site',help='Balsam Site on which to generate the integration grids.')
   parser.add_option('-i','--evgen-site',dest='evgen_site',help='Balsam Site on which to generate the events.')
   parser.add_option('-n','--minEvents',dest='evgen_site',help='Balsam Site on which to generate the events.')
   options,args = parser.parse_args()


   manditory_args = [
                     'ecmEnergy',
                     'runNumber',
                     'jobConfig',
                     'evgenJobOpts',
                     'outputEVNTFile',
                     'AMITag',
                     'steering',
                     'gridgen_site',
                     'evgen_site',
                  ]

   for man in manditory_args:
      if options.__dict__[man] is None:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)

   submit_alpgen_pilot(
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
             )




def submit_alpgen_pilot(
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
                 increment_iseed1      = True,
                 ):

   taskID=str(int(time.time()*1000000))
   
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
   if iseed2 is not None:
      job.input_iseed2 = iseed2
   if iseed3 is not None:
      job.input_iseed3 = iseed3
   if iseed4 is not None:
      job.input_iseed4 = iseed4
   
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
