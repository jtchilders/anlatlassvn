import os,sys,logging,shutil
logger = logging.getLogger(__name__)
from sherpa.SherpaInputCard import SherpaInputCard
import GridFtp
sys.path.append('/users/hpcusers/argobalsam/production/argobalsam/common_core')
from ArgoJob import ArgoJob
from BalsamJob import BalsamJob

class SherpaArgoJob:
   INPUT_FILE_CODE_GEN     = 'codegen.Run.dat'
   INPUT_FILE_INTEGRATION  = 'integration.Run.dat'
   INPUT_FILE_EVTGEN       = 'evtgen.Run.dat'
   PYTHON_INPUT_FILE       = 'Run.py'

   SHERPA_EVENT_OUTPUT     = 'EVENT_OUTPUT=HepMc_GenEvent'
   SHERPA_N_EVENTS         = 'EVENTS'

   CODEGEN_OPTS            = 'INIT_ONLY=1'
   INTEGRATION_OPTS        = 'AMEGIC_LIBRARY_MODE=1'
   INTEGRATION_BGQ_OPTS    = 'SQLITE_OPEN_FLAG=unix-none'
   EVTGEN_OPTS             = 'AMEGIC_LIBRARY_MODE=1'
   EVTGEN_BGQ_OPTS         = 'SQLITE_OPEN_FLAG=unix-none'


   
   CODEGEN_OUTPUT_TARBALL  = 'sherpa_codegen_output.tar.gz'
   CODEGEN_POSTSUBMIT      = 'sherpa_codegen_postsubmit.sh'
   CODEGEN_POSTSUBMIT_LOG  = 'sherpa_codegen_postsubmit.log'
   CODEGEN_LOGFILENAME     = 'codegen.log'
   
   INTEGRATION_OUTPUT_TARBALL  = 'sherpa_integration_output.tar.gz'
   INTEGRATION_PRESUBMIT   = 'sherpa_integration_presubmit.sh'
   INTEGRATION_PRESUBMIT_LOG = 'sherpa_integration_presubmit.sh.log'
   INTEGRATION_POSTSUBMIT  = 'sherpa_integration_postsubmit.sh'
   INTEGRATION_POSTSUBMIT_LOG = 'sherpa_integration_postsubmit.sh.log'
   INTEGRATION_LOGFILENAME = 'integration.log'
   
   EVTGEN_OUTPUT_TARBALL   = 'sherpa_evtgen_output.tar.gz'
   EVTGEN_PRESUBMIT        = 'sherpa_evtgen_presubmit.sh'
   EVTGEN_PRESUBMIT_LOG    = 'sherpa_evtgen_presubmit.sh.log'
   EVTGEN_POSTSUBMIT       = 'sherpa_evtgen_postsubmit.sh'
   EVTGEN_POSTSUBMIT_LOG   = 'sherpa_evtgen_postsubmit.sh.log'
   EVTGEN_LOGFILENAME      = 'evtgen.log'
   EVTGEN_HEPMC_BASE       = 'events'


   

   SHERPA_EXECUTABLE_SCRIPT = 'sherpa.sh'
   SHERPA_EXECUTABLE       = 'Sherpa'

   def __init__(self,
                code_gen_site                   = None,
                code_gen_input_filename         = None,

                integration_site                = None,
                integration_input_filename      = None,
                integration_nodes               = None,
                integration_processes_per_node  = None,
                integration_wall_minutes        = None,

                evtgen_site                     = None,
                evtgen_input_filename           = None,
                evtgen_nodes                    = None,
                evtgen_processes_per_node       = None,
                evtgen_wall_minutes             = None,
                evtgen_number_events            = None,

                working_path                    = None,
                input_url                       = None,
                output_url                      = None,
                username                        = None,
                random_seed1                    = None,
                random_seed2                    = None,
                random_seed3                    = None,
                random_seed4                    = None,
                python_input_card               = None,
                group_identifier                = None,
               ):
      self.code_gen_site                     = code_gen_site
      self.code_gen_input_filename           = code_gen_input_filename

      self.integration_site                  = integration_site
      self.integration_input_filename        = integration_input_filename
      self.integration_nodes                 = integration_nodes
      self.integration_processes_per_node    = integration_processes_per_node
      self.integration_wall_minutes          = integration_wall_minutes

      self.evtgen_site                       = evtgen_site
      self.evtgen_input_filename             = evtgen_input_filename
      self.evtgen_nodes                      = evtgen_nodes
      self.evtgen_processes_per_node         = evtgen_processes_per_node
      self.evtgen_wall_minutes               = evtgen_wall_minutes
      self.evtgen_number_events              = evtgen_number_events

      self.working_path                      = working_path
      self.input_url                         = input_url
      self.output_url                        = output_url
      self.username                          = username
      self.random_seed1                      = random_seed1
      self.random_seed2                      = random_seed2
      self.random_seed3                      = random_seed3
      self.random_seed4                      = random_seed4
      self.python_input_card                 = python_input_card
      self.group_identifier                  = group_identifier

   def get_argo_job(self):
      ##-----------------------
      # setup input files
      ##-----------------------

      if (self.code_gen_input_filename is None and
          self.integration_input_filename is None and
          self.evtgen_input_filename is None and
          self.python_input_card is None):
         logger.error('no input cards specified, nothing to work with.')
         return None

      if (self.code_gen_input_filename is not None and
             not os.path.exists(os.path.join(self.working_path,self.INPUT_FILE_CODE_GEN))):
         shutil.copyfile(self.code_gen_input_filename,os.path.join(self.working_path,self.INPUT_FILE_CODE_GEN))
      if (self.integration_input_filename is not None and
             not os.path.exists(os.path.join(self.working_path,self.INPUT_FILE_INTEGRATION))):
         shutil.copyfile(self.integration_input_filename,os.path.join(self.working_path,self.INPUT_FILE_INTEGRATION))
      if (self.evtgen_input_filename is not None and
             not os.path.exists(os.path.join(self.working_path,self.INPUT_FILE_EVTGEN))):
         shutil.copyfile(self.evtgen_input_filename,os.path.join(self.working_path,self.INPUT_FILE_EVTGEN))

      
      # # setup integration input card
      # if self.integration_input_filename is not None:
      #    integration_card = SherpaInputCard.read_file(self.integration_input_filename)
      #    if integration_card is None:
      #       logger.error(" ERROR reading input file for integration step: " + self.integration_input_filename)
      #       return None
      #    # update random seeds
      #    if ( self.random_seed1 is not None and 
      #         self.random_seed2 is not None and 
      #         self.random_seed3 is not None and 
      #         self.random_seed4 is not None
      #       ):
      #       integration_card.set_random_seeds(self.random_seed1,self.random_seed2,self.random_seed3,self.random_seed4)
      #    integration_card.write_file(os.path.join(self.working_path,self.INPUT_FILE_INTEGRATION))
      
      # # setup event generation input card
      # if self.evtgen_input_filename is not None:
      #    evtgen_card = SherpaInputCard.read_file(self.evtgen_input_filename)
      #    if evtgen_card is None:
      #       logger.error(" ERROR reading input file for event generation step: " + self.evtgen_input_filename)
      #       return None
      #    if self.evtgen_number_events is not None:
      #       evtgen_card.set_number_events(self.evtgen_number_events)
      #    # update random seeds
      #    if ( self.random_seed1 is not None and 
      #         self.random_seed2 is not None and 
      #         self.random_seed3 is not None and 
      #         self.random_seed4 is not None
      #       ):
      #       evtgen_card.set_random_seeds(self.random_seed1,self.random_seed2,self.random_seed3,self.random_seed4)
      #    evtgen_card.write_file(os.path.join(self.working_path,self.INPUT_FILE_EVTGEN))
      

      # common sherpa command line options
      common_opts = 'MPI_SEED_MODE=1'
      integration_opts = self.INTEGRATION_OPTS
      codegen_opts = self.CODEGEN_OPTS
      evtgen_opts = (self.EVTGEN_OPTS + ' ' 
                  + self.SHERPA_EVENT_OUTPUT + '[' + self.EVTGEN_HEPMC_BASE + ']')
      if self.evtgen_number_events is not None:
            evtgen_opts += ' ' + self.SHERPA_N_EVENTS + '=' + str(self.evtgen_number_events)

      if self.random_seed1 is not None:
         common_opts += ' RANDOM_SEED1=' + str(self.random_seed1)
      if self.random_seed2 is not None:
         common_opts += ' RANDOM_SEED2=' + str(self.random_seed2)
      if self.random_seed3 is not None:
         common_opts += ' RANDOM_SEED3=' + str(self.random_seed3)
      if self.random_seed4 is not None:
         common_opts += ' RANDOM_SEED4=' + str(self.random_seed4)

      if self.python_input_card is not None:
         common_opts += ' RUNDATA=' + self.PYTHON_INPUT_FILE

         common_opts += ' BEAM_1=2212 BEAM_2=2212 BEAM_ENERGY_1=6500 BEAM_ENERGY_2=6500 MASS[6]=172.5 MASS[23]=91.1876 MASS[24]=80.399 WIDTH[23]=2.4952 WIDTH[24]=2.085 WIDTH[15]=2.26735e-12 SIN2THETAW=0.23113 EVENT_GENERATION_MODE=Weighted FRAGMENTATION=Off MI_HANDLER=None'
         codegen_opts += ' LOG_FILE=' + self.CODEGEN_LOGFILENAME
         integration_opts += ' EVENTS=1 LOG_FILE=' + self.INTEGRATION_LOGFILENAME
         evtgen_opts += ' LOG_FILE=' + self.EVTGEN_LOGFILENAME
      else:
         codegen_opts += ' -l ' + self.CODEGEN_LOGFILENAME
         integration_opts += ' -l ' + self.INTEGRATION_LOGFILENAME
         evtgen_opts += ' -l ' + self.EVTGEN_LOGFILENAME




      # create code generation balsam job
      code_gen = None
      if self.code_gen_site is not None:
         # create input card if needed
         if self.python_input_card is not None:
            if os.path.exists(self.python_input_card):
               shutil.copyfile(self.python_input_card,os.path.join(self.working_path,self.PYTHON_INPUT_FILE))
            else:
               logger.error(' File not found: ' + self.python_input_card)
               return None
         else:
            code_gen_card = SherpaInputCard.read_file(self.code_gen_input_filename)
            if code_gen_card is None:
               logger.error(" ERROR reading input file for code generation step: " + self.code_gen_input_filename)
               return None
            code_gen_card.write_file(os.path.join(self.working_path,self.INPUT_FILE_CODE_GEN))

         code_gen = BalsamJob()
         
         code_gen.target_site                = self.code_gen_site
         
         if code_gen.target_site == 'argo_cluster':
            code_gen.executable                 = self.SHERPA_EXECUTABLE_SCRIPT
         else:
            code_gen.executable                 = self.SHERPA_EXECUTABLE
         
         code_gen.executable_args               = codegen_opts + ' ' + common_opts
         if self.python_input_card is None:
            code_gen.executable_args            = '-f ' + self.INPUT_FILE_CODE_GEN + ' ' + code_gen.executable_args
         
         code_gen.input_files                   = []
         if self.python_input_card is not None:
            code_gen.input_files += [self.PYTHON_INPUT_FILE]
         else:
            code_gen.input_files += [self.INPUT_FILE_CODE_GEN]
         
         code_gen.output_files                  = [SherpaArgoJob.CODEGEN_OUTPUT_TARBALL,
                                                   SherpaArgoJob.CODEGEN_POSTSUBMIT_LOG,
                                                   self.CODEGEN_LOGFILENAME,
                                                  ]
         
         code_gen.nodes                         = 1 # code generation is serial, fails with more than one rank
         code_gen.processes_per_node            = 1 # code generation is serial, fails with more than one rank
         if code_gen.target_site == 'argo_cluster':
            code_gen.executable_args = str(code_gen.processes_per_node) + ' ' + code_gen.executable_args # need number of nodes, in this case 1
         
         code_gen.wall_minutes                  = 0 # serial, so runs on argo_cluster
         code_gen.username                      = self.username
         code_gen.postprocess                   = SherpaArgoJob.CODEGEN_POSTSUBMIT
         code_gen.postprocess_args              = SherpaArgoJob.CODEGEN_OUTPUT_TARBALL
      
      # create integration balsam job
      integration = None
      if self.integration_site is not None:
         integration = BalsamJob()
         
         integration.target_site              = self.integration_site
         
         integration.executable_args          = ' ' + integration_opts + ' ' + common_opts

         # options that try to optimize the integration
         # integration should do 10 events per rank per interation
         integration.executable_args          = ' PSI_ITMIN='+str(self.integration_nodes*self.integration_processes_per_node*10)
         # do not perform 30 iterations if it is not needed
         integration.executable_args          = ' FINISH_OPTIMIZATION=Off'

         if integration.target_site in ['argo_cluster','edison','mira','cetus','vesta']:
            integration.executable            = self.SHERPA_EXECUTABLE_SCRIPT
         else:
            integration.executable            = self.SHERPA_EXECUTABLE
         
         if self.python_input_card is None:
            integration.executable_args            = '-f ' + self.INPUT_FILE_INTEGRATION + ' ' + integration.executable_args
         
         if integration.target_site in ['vesta','cetus','mira']:
            integration.executable_args       += ' ' + self.INTEGRATION_BGQ_OPTS
            integration.scheduler_args        = ' --mode=script'
         
         integration.input_files              = [SherpaArgoJob.CODEGEN_OUTPUT_TARBALL]
         if self.python_input_card is not None:
            integration.input_files += [self.PYTHON_INPUT_FILE]
         else:
            integration.input_files += [self.INPUT_FILE_INTEGRATION]
         
         integration.output_files             = [SherpaArgoJob.INTEGRATION_OUTPUT_TARBALL,
                                                 SherpaArgoJob.INTEGRATION_PRESUBMIT_LOG,
                                                 SherpaArgoJob.INTEGRATION_POSTSUBMIT_LOG,
                                                 self.INTEGRATION_LOGFILENAME,
                                                ]
         
         integration.nodes                    = self.integration_nodes
         integration.processes_per_node       = self.integration_processes_per_node
         if integration.target_site in [ 'argo_cluster','mira','cetus','vesta']:
            integration.executable_args       = ' ' + str(integration.processes_per_node) + ' ' + integration.executable_args
         

         integration.wall_minutes             = self.integration_wall_minutes
         integration.username                 = self.username
         integration.preprocess               = SherpaArgoJob.INTEGRATION_PRESUBMIT
         integration.preprocess_args          = SherpaArgoJob.CODEGEN_OUTPUT_TARBALL
         integration.postprocess              = SherpaArgoJob.INTEGRATION_POSTSUBMIT
         integration.postprocess_args         = SherpaArgoJob.INTEGRATION_OUTPUT_TARBALL

         if integration.target_site == 'edison':
            integration.scheduler_args = " -v RANKS_PER_NODE=" + str(integration.processes_per_node) + ",SHERPA_ARGS='" + integration.executable_args + "'"
            integration.executable_args = None


      # create event generation job
      evtgen = None
      if self.evtgen_site is not None:
         evtgen = BalsamJob()
         
         evtgen.target_site              = self.evtgen_site
         
         if evtgen.target_site in ['argo_cluster','edison','mira','cetus','vesta']:
            evtgen.executable            = self.SHERPA_EXECUTABLE_SCRIPT
         else:
            evtgen.executable            = self.SHERPA_EXECUTABLE
         
         evtgen.executable_args          = ' ' + evtgen_opts + ' ' + common_opts
         if self.python_input_card is None:
            evtgen.executable_args            = '-f ' + self.INPUT_FILE_EVTGEN + ' ' + evtgen.executable_args
         
         if evtgen.target_site in [ 'vesta','cetus','mira']:
            evtgen.executable_args       += ' ' + self.EVTGEN_BGQ_OPTS
            integration.scheduler_args        = ' --mode=script'
         
         evtgen.input_files              = [SherpaArgoJob.INTEGRATION_OUTPUT_TARBALL]
         if self.python_input_card is not None:
            evtgen.input_files += [self.PYTHON_INPUT_FILE]
         else:
            evtgen.input_files += [self.INPUT_FILE_EVTGEN]
         
         evtgen.output_files             = [self.EVTGEN_OUTPUT_TARBALL,
                                            self.EVTGEN_PRESUBMIT_LOG,
                                            self.EVTGEN_POSTSUBMIT_LOG,
                                            self.EVTGEN_LOGFILENAME,
                                           ]
         
         evtgen.nodes                    = self.evtgen_nodes
         evtgen.processes_per_node       = self.evtgen_processes_per_node
         if evtgen.target_site in ['argo_cluster','mira','cetus','vesta']:
            evtgen.executable_args       = ' ' + str(evtgen.processes_per_node) + ' ' + evtgen.executable_args
         

         evtgen.wall_minutes             = self.evtgen_wall_minutes
         evtgen.username                 = self.username
         evtgen.preprocess               = self.EVTGEN_PRESUBMIT
         evtgen.preprocess_args          = self.INTEGRATION_OUTPUT_TARBALL
         evtgen.postprocess              = self.EVTGEN_POSTSUBMIT
         evtgen.postprocess_args         = self.EVTGEN_OUTPUT_TARBALL

         if evtgen.target_site == 'edison':
            evtgen.scheduler_args = " -v RANKS_PER_NODE=" + str(evtgen.processes_per_node) + ",SHERPA_ARGS='" + evtgen.executable_args + "'"
            evtgen.executable_args = None
      
      argojob = ArgoJob()
      argojob.input_url          = self.input_url
      argojob.output_url         = self.output_url
      argojob.username           = self.username
      argojob.group_identifier   = self.group_identifier
      if code_gen is not None:
         argojob.add_job(code_gen)
      if integration is not None:
         argojob.add_job(integration)
      if evtgen is not None:
         argojob.add_job(evtgen)


      # copy files to grid ftp location
      try:
         GridFtp.globus_url_copy(self.working_path + '/',self.input_url + '/')
      except:
         logger.error(' received exception while copying working path to grid ftp input path: ' + str(sys.exc_info()[1]))
         raise

      return argojob

