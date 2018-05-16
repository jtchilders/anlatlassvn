import os,sys,logging,shutil
logger = logging.getLogger(__name__)
sys.path.append('/users/hpcusers/argobalsam/production/argobalsam/common_core')
from ArgoJob import ArgoJob
from BalsamJob import BalsamJob
import GridFtp

class PythiaArgoJob:
   
   PYTHIA_HEPMC_FILEBASE      = 'pythia_output'
   OUTPUT_TARBALL             = PYTHIA_HEPMC_FILEBASE + '.gz.tar'
   RANKS_PER_COLLECTOR        = 128
   EVTGEN_SCHEDULER_ARGS      = '--mode=script'
   EXECUTABLE                 = 'pythia'
   PREPROCESS                 = 'pythia_presubmit.sh'
   POSTPROCESS                = 'pythia_postsubmit.py'

   def __init__(self,
                input_filename                  = None,
                working_path                    = None,
                nevts                           = None,
                site                            = None,
                nodes                           = None,
                ranks_per_node                  = None,
                input_url                       = None,
                output_url                      = None,
                username                        = None,
                random_seed                     = None,
                wall_minutes                    = None,
                redirect_std                    = False,
                disable_output                  = False,
                gzip_output                     = True,
               ):
      self.input_filename                    = input_filename
      self.working_path                      = working_path
      self.nevts                             = nevts
      self.site                              = site
      self.nodes                             = nodes
      self.ranks_per_node                    = ranks_per_node
      self.input_url                         = input_url
      self.output_url                        = output_url
      self.username                          = username
      self.random_seed                       = random_seed
      self.wall_minutes                      = wall_minutes
      self.redirect_std                      = redirect_std
      self.disable_output                    = disable_output
      self.gzip_output                       = gzip_output

   def get_argo_job(self):
      
       
      # calculate the number of events per rank to calculate
      evts_per_rank = int(self.nevts/(self.nodes*self.ranks_per_node))
      if evts_per_rank < 10:
         logger.error('Too few events per rank = ' + str(evts_per_rank))
         return None

      logger.info('Each rank going to produce ' + str(evts_per_rank) + ' events.')
   
      # build pythia command line args
      exe_args  = ' -n ' + str(evts_per_rank)
      exe_args += ' -a ' + str(PythiaArgoJob.RANKS_PER_COLLECTOR)
      exe_args += ' -x xmldoc'
      exe_args += ' -c ' + os.path.basename(self.input_filename)
      exe_args += ' -o ' + PythiaArgoJob.PYTHIA_HEPMC_FILEBASE
      if self.gzip_output:
         exe_args += ' -g ' # gzip output
      if self.random_seed is not None:
         exe_args += ' -s ' + str(self.random_seed)
      if self.redirect_std:
         exe_args += ' -r ' # redirect stdout/stderr output to files
      if self.disable_output:
         exe_args += ' -d ' # disable data output to file

      ##-----------------------
      # setup input files
      ##-----------------------
      
      # copy input file to working path
      try:
         shutil.copy(self.input_filename,self.working_path + '/')
      except:
         logger.exception(' received exception while copying input file: ' + str(sys.exc_info()[1]))

      # copy files to grid ftp location
      try:
         GridFtp.globus_url_copy(self.working_path + '/',self.input_url + '/')
      except:
         logger.exception(' received exception while copying working path to grid ftp input path: ' + str(sys.exc_info()[1]))
      
      # create executable
      executable = PythiaArgoJob.EXECUTABLE

      # create warmup balsam job
      balsam_job = BalsamJob()
      balsam_job.executable         = executable
      balsam_job.executable_args    = exe_args
      balsam_job.input_files        = [self.input_filename]
      balsam_job.output_files       = [PythiaArgoJob.OUTPUT_TARBALL]
      balsam_job.target_site        = self.site
      balsam_job.nodes              = self.nodes
      balsam_job.processes_per_node = self.ranks_per_node
      balsam_job.wall_minutes       = self.wall_minutes
      balsam_job.username           = self.username
      balsam_job.preprocess         = PythiaArgoJob.PREPROCESS
      balsam_job.postprocess        = PythiaArgoJob.POSTPROCESS
      
            
      
      argojob = ArgoJob()
      argojob.input_url          = self.input_url
      argojob.output_url         = self.output_url
      argojob.username           = self.username
      argojob.add_job(balsam_job)

      return argojob

