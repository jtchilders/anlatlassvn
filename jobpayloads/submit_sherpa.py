#!/usr/bin/env python
import os,sys,logging,optparse,time,logging,json,shutil,glob

sys.path.append('/users/hpcusers/svn/generators/alpgen/v214/usercode/python')
try:

 from ArgoJobs.SherpaArgoJob import SherpaArgoJob
except:
  import traceback
  traceback.print_exc()
  sys.exit(1)
sys.path.append('/users/hpcusers/balsam_production/balsam_deploy/common_core')
from MessageInterface import MessageInterface
from ArgoJob import ArgoJob

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

DEFAULT_CODEGEN_SITE = 'argo_cluster'
JOB_DEPOT = 'jobs_sherpa'

def main():
   parser = optparse.OptionParser(description='Submit Sherpa job to ARGO')

   # code generation step options
   parser.add_option('-a','--codegen-card',dest='codegen_card',help='Code Generation Step: Sherpa configuration card for code generation step. If not specified, code generation step is skipped.')
   parser.add_option('-b','--codegen-site',dest='codegen_site',help='Code Generation Step: Site on which to run. For the moment, this step is serial only. ')

   # integration step options
   parser.add_option('-c','--int-card',dest='integration_card',help='Integration Step: The Sherpa input file which carries all the options. If not specified, integration step is skipped.')
   parser.add_option('-e','--int-site',dest='integration_site',help='Integration Step: Site on which to run.')
   parser.add_option('-f','--int-num-nodes',dest='integration_numnodes',help='Integration Step: Number of nodes to use on destination site',type='int')
   parser.add_option('-g','--int-cpu-per-node',dest='integration_cpu_per_node',help='Integration Step: Number of CPUs per node to use on destination site',type='int')
   parser.add_option('-i','--int-walltime',dest='integration_walltime',help='Integration Step: The wall time to use when submitting to the queue (in minutes).',type='int')

   # event generation step options
   parser.add_option('-j','--evtgen-card',dest='evtgen_card',help='Event Generation Step: The Sherpa input file which carries all the options. If not specified, integration step is skipped.')
   parser.add_option('-k','--evtgen-site',dest='evtgen_site',help='Event Generation Step: Site on which to run.')
   parser.add_option('-l','--evtgen-num-nodes',dest='evtgen_numnodes',help='Event Generation Step: Number of nodes to use on destination site',type='int')
   parser.add_option('-m','--evtgen-cpu-per-node',dest='evtgen_cpu_per_node',help='Event Generation Step: Number of CPUs per node to use on destination site',type='int')
   parser.add_option('-n','--evtgen-walltime',dest='evtgen_walltime',help='Event Generation Step: The wall time to use when submitting to the queue (in minutes).',type='int')
   parser.add_option('-o','--evtgen-nevts',dest='evtgen_nevts',help='Event Generation Step: The number of events to generate per thread of Sherpa.',type='int')

   parser.add_option('-p','--group-identifier',dest='group_identifier',help='An identifier to label groups of jobs.')
   
   parser.add_option('-r','--regen',dest='regen_path',help='Regenerate from an earlier job. If this path is given and the code generation input file is not given, but the input file for the other two steps are given, this path will be used as the code source for the integration and event generation. If the code generation and integration input files are given, then this path will be used for the event generation step.')

   parser.add_option('-x','--no-submit',dest='submit',help='do not submit the message to ARGO. For testing purposes.',action='store_false',default=True)

   parser.add_option('','--randomseed1',dest='randomseed1',help='override random seed 1 in sherpa input file. For the moment, all four must be provided together.',type='int')
   parser.add_option('','--randomseed2',dest='randomseed2',help='override random seed 2 in sherpa input file. For the moment, all four must be provided together.',type='int')
   parser.add_option('','--randomseed3',dest='randomseed3',help='override random seed 3 in sherpa input file. For the moment, all four must be provided together.',type='int')
   parser.add_option('','--randomseed4',dest='randomseed4',help='override random seed 4 in sherpa input file. For the moment, all four must be provided together.',type='int')

   parser.add_option('','--python-input-card',dest='python_input_card',help='Sherpa can extract run card from python file. Using this option overrides using the codegen-card, integration-card, and evtgen-card settings and makes Sherpa use this python file.')

   options,args = parser.parse_args()
   
   submit_sherpa(
                 options.codegen_card,
                 options.codegen_site,

                 options.integration_card,
                 options.integration_site,
                 options.integration_numnodes,
                 options.integration_cpu_per_node,
                 options.integration_walltime,

                 options.evtgen_card,
                 options.evtgen_site,
                 options.evtgen_numnodes,
                 options.evtgen_cpu_per_node,
                 options.evtgen_walltime,
                 options.evtgen_nevts,

                 options.regen_path,
                 options.submit,
                 options.group_identifier,

                 options.randomseed1,
                 options.randomseed2,
                 options.randomseed3,
                 options.randomseed4,

                 options.python_input_card,
                )




def submit_sherpa(
                 codegen_card,
                 codegen_site,

                 integration_card,
                 integration_site,
                 integration_numnodes,
                 integration_cpu_per_node,
                 integration_wall_minutes,

                 evtgen_card,
                 evtgen_site,
                 evtgen_numnodes,
                 evtgen_cpu_per_node,
                 evtgen_wall_minutes,
                 evtgen_number_events,

                 regen_path                  = None,
                 submit                      = True,
                 group_identifier            = None,

                 randomseed1                 = None,
                 randomseed2                 = None,
                 randomseed3                 = None,
                 randomseed4                 = None,

                 python_input_card           = None,
                 ):

   # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   #  Input Consistency Checks
   # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   

   # ensure all 4 random seeds are given or not given
   if (    not (randomseed1 is None and 
                randomseed2 is None and 
                randomseed3 is None and 
                randomseed4 is None) 
       and not (randomseed1 is not None and 
                randomseed2 is not None and 
                randomseed3 is not None and 
                randomseed4 is not None)
      ):
      logger.error('For the moment, all random seeds must be specified together or not at all.')
      sys.exit(-1)

   # sanity checks on Mira partition sizes
   if integration_site == 'mira':
      numnode = integration_numnodes
      if( numnode != 2**9
          and numnode != 2**10
          and numnodes != 2**11
          and numnodes != 2**12
          and numnodes != 2**13
          and numnodes != 2**14
          and numnodes != 2**15
          and numnodes != (2**15 + 2**14)
        ):
        parser.error(' You selected Mira as a site and are using a number of nodes, ' + str(numnodes) + ', that is not a power of 2 (512,1024,2048,4096,8192...)')
        sys.exit(-1)
   if evtgen_site == 'mira':
      numnode = evtgen_numnodes
      if( numnode != 2**9
          and numnode != 2**10
          and numnodes != 2**11
          and numnodes != 2**12
          and numnodes != 2**13
          and numnodes != 2**14
          and numnodes != 2**15
          and numnodes != (2**15 + 2**14)
        ):
        parser.error(' You selected Mira as a site and are using a number of nodes, ' + str(numnodes) + ', that is not a power of 2 (512,1024,2048,4096,8192...)')
        sys.exit(-1)

   if regen_path is not None:
      if not os.path.exists(regen_path):
         logger.error('Regeneration path given, but does not exist: ' + regen_path)
         sys.exit(-1)


   # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   #  Sherpa Job Preparation
   # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

   TOP_PATH = os.getcwd() # directory in which script was run
   
   job_id = create_job_path()
   logger.info('job_id: ' + str(job_id))

   RUN_PATH = os.path.join(TOP_PATH,JOB_DEPOT,str(job_id)) # directory in which to store files
   logger.info('run path: ' + RUN_PATH)

   if codegen_card is not None:
      shutil.copyfile(codegen_card,os.path.join(RUN_PATH,os.path.basename(codegen_card)))
   if integration_card is not None:
      shutil.copyfile(integration_card,os.path.join(RUN_PATH,os.path.basename(integration_card)))
   if evtgen_card is not None:
      shutil.copyfile(evtgen_card,os.path.join(RUN_PATH,os.path.basename(evtgen_card)))
   if python_input_card is not None:
      shutil.copyfile(python_input_card,os.path.join(RUN_PATH,os.path.basename(python_input_card)))
   
   if regen_path is not None:
      logger.info('regenerating from path: ' + regen_path)
      # must first delete the directory we just created
      #os.rmdir(RUN_PATH)
      # because copytree must create the destination directory
      #shutil.copytree(regen_path,RUN_PATH)
      possible_files_to_copy = []
      if codegen_site is not None:
         logger.info(" regenerate path specificed and codegen specified. That doesn't make sense.... exiting.")
         sys.exit(-1)
      
      possible_files_to_copy.append(SherpaArgoJob.INPUT_FILE_CODE_GEN)
      possible_files_to_copy.append(SherpaArgoJob.CODEGEN_OUTPUT_TARBALL)
      possible_files_to_copy.append(SherpaArgoJob.CODEGEN_LOGFILENAME)
      possible_files_to_copy.append(SherpaArgoJob.CODEGEN_POSTSUBMIT_LOG)
      possible_files_to_copy.append(SherpaArgoJob.PYTHON_INPUT_FILE)

      if integration_site is None:
         possible_files_to_copy.append(SherpaArgoJob.INPUT_FILE_INTEGRATION)
         possible_files_to_copy.append(SherpaArgoJob.INTEGRATION_OUTPUT_TARBALL)
         possible_files_to_copy.append(SherpaArgoJob.INTEGRATION_LOGFILENAME)
         possible_files_to_copy.append(SherpaArgoJob.INTEGRATION_PRESUBMIT_LOG)
         possible_files_to_copy.append(SherpaArgoJob.INTEGRATION_POSTSUBMIT_LOG)

      files_to_copy = [ file for file in possible_files_to_copy if os.path.exists(os.path.join(regen_path,file)) ]
      for file in files_to_copy:
         shutil.copyfile(os.path.join(regen_path,file),os.path.join(RUN_PATH,file))

   cwd = os.getcwd()
   os.chdir(RUN_PATH)

   # create input/output gridftp urls
   grid_ftp_server = 'atlasgridftp02.hep.anl.gov'
   grid_base_path  = '/grid/atlas/hpc/argo/jobs/' + str(job_id)
   input_url = 'gsiftp://' + grid_ftp_server + grid_base_path
   output_url = input_url

   job = SherpaArgoJob(
               code_gen_site                   = codegen_site,
               code_gen_input_filename         = codegen_card,

               integration_site                = integration_site,
               integration_input_filename      = integration_card,
               integration_nodes               = integration_numnodes,
               integration_processes_per_node  = integration_cpu_per_node,
               integration_wall_minutes        = integration_wall_minutes,

               evtgen_site                     = evtgen_site,
               evtgen_input_filename           = evtgen_card,
               evtgen_nodes                    = evtgen_numnodes,
               evtgen_processes_per_node       = evtgen_cpu_per_node,
               evtgen_wall_minutes             = evtgen_wall_minutes,
               evtgen_number_events            = evtgen_number_events,

               working_path                    = RUN_PATH,
               input_url                       = input_url,
               output_url                      = output_url,
               username                        = os.environ['USER'],
               random_seed1                    = randomseed1,
               random_seed2                    = randomseed2,
               random_seed3                    = randomseed3,
               random_seed4                    = randomseed4,
               python_input_card               = python_input_card,
               group_identifier                = group_identifier,
            )
   
   argojob = job.get_argo_job()
   if argojob is None:
      logger.error(' error getting argo job ')
      sys.exit(-10)
   #argojob.email_address = 'jchilders@anl.gov,turam@anl.gov'
   #argojob.email_address = 'turam@anl.gov'
   argojob.email_address = 'jchilders@anl.gov'
   
   job_txt = argojob.serialize()


   #print "JOB TEXT"
   #print job_txt

   if submit:
      mi = MessageInterface()
      mi.host = 'atlasgridftp02.hep.anl.gov'
      mi.port = 5671
      mi.ssl_cert = os.environ['X509_USER_CERT'] #'/users/hpcusers/balsam/gridsecurity/jchilders/xrootdsrv-cert.pem'
      mi.ssl_key  = os.environ['X509_USER_KEY'] #'/users/hpcusers/balsam/gridsecurity/jchilders/xrootdsrv-key.pem'
      mi.ssl_ca_certs = os.environ['X509_CA_CERTS'] #'/users/hpcusers/balsam/gridsecurity/jchilders/cacerts.pem'
      mi.exchange_name = 'argo_users'
      #if dev:
      #   mi.exchange_name = 'argo_users_dev'
      print 'opening connection'
      mi.open_blocking_connection()
      routing_key = 'argo_job'
      #if dev:
      #   routing_key = 'argo_job_dev'
      print ' sending msg '
      mi.send_msg(argojob.serialize(),routing_key)
      print ' done sending'
      mi.close()
      print ' closing connection'
   else:
      logger.info(' not submitting job\n'+job_txt)

   os.chdir(cwd)

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
