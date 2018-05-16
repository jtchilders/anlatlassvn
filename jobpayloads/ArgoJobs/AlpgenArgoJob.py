import os,sys,logging,shutil
logger = logging.getLogger(__name__)
from AlpgenInputFile import AlpgenInputFile
import GridFtp
sys.path.append('/users/hpcusers/argobalsam/production/argobalsam/common_core')
from ArgoJob import ArgoJob
from BalsamJob import BalsamJob

class AlpgenArgoJob:
   INPUT_FILE_POSTFIX_IMODE0 = '.input.0'
   INPUT_FILE_POSTFIX_IMODE1 = '.input.1'
   INPUT_FILE_POSTFIX_IMODE2 = '.input.2'

   EVTGEN_EXECUTABLE         = 'alpgenCombo.sh'

   def __init__(self,
                process                         = None,
                input_filename                  = None,
                input_iseed1                    = None,
                input_iseed2                    = None,
                input_iseed3                    = None,
                input_iseed4                    = None,
                warmup_phase0_number_events     = None,
                warmup_phase0_number_iterations = None,
                warmup_phase1_number_events     = None,
                warmup_wall_minutes             = None,
                evtgen_phase0_number_events     = None,
                evtgen_phase0_number_iterations = None,
                evtgen_phase1_number_events     = None,
                evtgen_nodes                    = None,
                evtgen_processes_per_node       = None,
                evtgen_wall_minutes             = None,
                working_path                    = None,
                input_url                       = None,
                output_url                      = None,
                pdf_filename                    = None,
                username                        = None,
                site                            = None,
                group_identifier                = None,
               ):
      self.process                           = process
      self.input_filename                    = input_filename
      self.input_iseed1                      = input_iseed1
      self.input_iseed2                      = input_iseed2
      self.input_iseed3                      = input_iseed3
      self.input_iseed4                      = input_iseed4
      self.warmup_phase0_number_events       = warmup_phase0_number_events
      self.warmup_phase0_number_iterations   = warmup_phase0_number_iterations
      self.warmup_phase1_number_events       = warmup_phase1_number_events
      self.warmup_wall_minutes               = warmup_wall_minutes
      self.evtgen_phase0_number_events       = evtgen_phase0_number_events
      self.evtgen_phase0_number_iterations   = evtgen_phase0_number_iterations
      self.evtgen_phase1_number_events       = evtgen_phase1_number_events
      self.evtgen_nodes                      = evtgen_nodes
      self.evtgen_processes_per_node         = evtgen_processes_per_node
      self.evtgen_wall_minutes               = evtgen_wall_minutes
      self.working_path                      = working_path
      self.input_url                         = input_url
      self.output_url                        = output_url
      self.pdf_filename                      = pdf_filename
      self.username                          = username
      self.warmup_site                       = 'argo_cluster'
      self.evtgen_site                       = site
      self.group_identifier                  = group_identifier

   def get_argo_job(self):
      ##-----------------------
      # setup input files
      ##-----------------------
      
      # load input file
      input = AlpgenInputFile()
      if input.read(self.input_filename) != 0:
         logger.error(" ERROR reading input file: " + self.input_filename)
         return None

      filename_base = input.filename_base
      
      # update random seeds
      if self.input_iseed1 is not None:
         logger.info('setting iseed1 = ' + str(self.input_iseed1))
         input.edit_option('iseed1',self.input_iseed1)
      if self.input_iseed2 is not None:
         logger.info('setting iseed2 = ' + str(self.input_iseed2))
         input.edit_option('iseed2',self.input_iseed2)
      if self.input_iseed3 is not None:
         logger.info('setting iseed3 = ' + str(self.input_iseed3))
         input.edit_option('iseed3',self.input_iseed3)
      if self.input_iseed4 is not None:
         logger.info('setting iseed4 = ' + str(self.input_iseed4))
         input.edit_option('iseed4',self.input_iseed4)

      # create input for imode 0
      input.imode             = 0
      input.start_with        = 0
      input.nevt              = self.warmup_phase0_number_events
      input.nitr              = self.warmup_phase0_number_iterations
      input.last_nevt         = self.warmup_phase1_number_events
      input_filename_imode0   = filename_base + self.INPUT_FILE_POSTFIX_IMODE0
      input.write(os.path.join(self.working_path,input_filename_imode0))
      
      # create input for imode 1
      input.imode             = 1
      input.start_with        = 2
      input.nevt              = self.evtgen_phase0_number_events
      input.nitr              = self.evtgen_phase0_number_iterations
      input.last_nevt         = self.evtgen_phase1_number_events
      input_filename_imode1   = filename_base + self.INPUT_FILE_POSTFIX_IMODE1
      input.write(os.path.join(self.working_path,input_filename_imode1))
      
      # create input for imode 2
      input.imode             = 2
      input.start_with        = 1
      input.nevt              = 0
      input.nitr              = 0
      input.last_nevt         = 0
      input_filename_imode2   = filename_base + self.INPUT_FILE_POSTFIX_IMODE2
      input.write(os.path.join(self.working_path,input_filename_imode2))
      
      # copy pdf file to working path
      try:
         shutil.copy(self.pdf_filename,self.working_path + '/')
      except:
         logger.error(' received exception while copying PDF file: ' + str(sys.exc_info()[1]))
         raise

      # copy files to grid ftp location
      try:
         logger.info(' GridFTP: ' + str(GridFtp.__file__))
         GridFtp.globus_url_copy(self.working_path + '/',self.input_url + '/')
      except:
         logger.error(' received exception while copying working path to grid ftp input path: ' + str(sys.exc_info()[1]))
         raise
      
      # create grid filenames
      grid1 = filename_base + '.grid1'
      grid2 = filename_base + '.grid2'

      # create executable
      executable = self.process + 'gen90_mpi'
      if any(x in self.evtgen_site for x in ['mira','vesta','cetus']):
         executable = self.process + 'gen90_mpi_nomrstpdfs'
      if 'edison' in self.evtgen_site:
         executable = self.process + 'gen90_mpi_nomrstpdfs'
      
      # create warmup balsam job
      warmup = BalsamJob()
      warmup.executable          = self.process + 'gen90_mpi'
      warmup.executable_args     = input_filename_imode0
      warmup.input_files         = [input_filename_imode0,
                                    os.path.basename(self.pdf_filename)]
      warmup.output_files        = [grid1,grid2]
      warmup.nodes               = 1
      warmup.processes_per_node  = 1
      warmup.wall_minutes        = self.warmup_wall_minutes
      warmup.username            = self.username
      warmup.target_site         = self.warmup_site
      
      # create filenames
      unw      = filename_base + '.unw'
      unw_gz   = unw + '.gz'
      unw_par  = filename_base + '_unw.par'
      wgt      = filename_base + '.wgt'
      wgt_par  = filename_base + '.par'
      directoryList_before = 'directoryList_before.txt'
      directoryList_after  = 'directoryList_after.txt'

      preprocess = 'alpgen_presubmit.sh'
      postprocess = 'alpgen_postsubmit.sh'

      # create event gen balsam job
      evtgen = BalsamJob()
      evtgen.executable          = self.EVTGEN_EXECUTABLE
      evtgen.executable_args     = executable + ' ' + input_filename_imode1 + ' ' + input_filename_imode2 + ' ' + str(self.evtgen_processes_per_node)
      evtgen.input_files         = [grid1,
                                    grid2,
                                    input_filename_imode1,
                                    input_filename_imode2,
                                    os.path.basename(self.pdf_filename)]
      evtgen.output_files        = [unw_gz,
                                    unw_par,
                                    #wgt,
                                    #wgt_par,
                                    directoryList_before,
                                    directoryList_after,
                                    postprocess + '.out',
                                    postprocess + '.err',
                                   ]
      evtgen.preprocess          = preprocess
      evtgen.postprocess         = postprocess
      evtgen.postprocess_args    = filename_base
      evtgen.nodes               = self.evtgen_nodes
      evtgen.processes_per_node  = self.evtgen_processes_per_node
      if 'argo_cluster' in self.evtgen_site:
         evtgen.processes_per_node = 1
      evtgen.wall_minutes        = self.evtgen_wall_minutes
      evtgen.username            = self.username
      if any(x in self.evtgen_site for x in ['mira','vesta','cetus']):
         evtgen.scheduler_args      = '--mode=script'
      #elif 'edison' in self.evtgen_site:
      #   evtgen.scheduler_args      = ('-v INPUT1=' + input_filename_imode1 + 
      #                                 ',INPUT2=' + input_filename_imode2 +
      #                                 ',EXE=' + executable +
      #                                 ',RANKS_PER_NODE=' + str(self.evtgen_processes_per_node)
      #                                )
      evtgen.target_site         = self.evtgen_site
      
      
      
      argojob = ArgoJob()
      argojob.input_url          = self.input_url
      argojob.output_url         = self.output_url
      argojob.username           = self.username
      argojob.group_identifier   = self.group_identifier
      argojob.add_job(warmup)
      argojob.add_job(evtgen)

      return argojob

