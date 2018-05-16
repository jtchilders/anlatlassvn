import os,sys,logging,shutil
logger = logging.getLogger(__name__)
sys.path.append('/users/hpcusers/balsam/balsam_deploy/common_core')
from ArgoJob import ArgoJob
from BalsamJob import BalsamJob
import GridFtp

class CondorArgoJob:
   def __init__(self,
                condor_job_file                 = None,
                condor_dagman_file              = None,
                working_path                    = None,
                #site                            = None,
                #nodes                           = None,
                #ranks_per_node                  = None,
                input_files                     = [],
                output_files                    = None,
                input_url                       = [],
                output_url                      = None,
                username                        = None,
                #wall_minutes                    = None,
               ):
      self.condor_job_file                   = condor_job_file
      self.condor_dagman_file                = condor_dagman_file
      self.working_path                      = working_path
      self.site                              = 'argo_cluster'
      self.nodes                             = 1
      self.ranks_per_node                    = 1
      self.input_files                       = input_files
      self.output_files                      = output_files
      self.input_url                         = input_url
      self.output_url                        = output_url
      self.username                          = username
      self.wall_minutes                      = -1

   def get_argo_job(self):
     
      ##-----------------------
      # setup input files
      ##-----------------------
      
      # copy input file to working path
      try:
         if self.condor_job_file:
            shutil.copy(self.condor_job_file,self.working_path + '/')
         if self.condor_dagman_file:
            shutil.copy(self.condor_dagman_file,self.working_path + '/')
         for file in self.input_files:
            shutil.copy(file,self.working_path + '/')
      except:
         logger.exception(' received exception while copying condor files ')
         raise

      # copy files to grid ftp location
      try:
         GridFtp.globus_url_copy(self.working_path + '/',self.input_url + '/')
      except:
         logger.exception(' received exception while copying working path to grid ftp input path: ' + str(sys.exc_info()[1]))
           # create warmup balsam job

      balsam_job = BalsamJob()
      balsam_job.condor_job_file    = self.condor_job_file
      balsam_job.condor_dagman_file = self.condor_dagman_file
      balsam_job.input_files        = self.input_files
      balsam_job.output_files       = self.output_files
      balsam_job.target_site        = self.site
      balsam_job.nodes              = self.nodes
      balsam_job.processes_per_node = self.ranks_per_node
      balsam_job.wall_minutes       = self.wall_minutes
      balsam_job.username           = self.username
      
      argojob = ArgoJob()
      argojob.input_url          = self.input_url
      argojob.output_url         = self.output_url
      argojob.username           = self.username
      argojob.add_job(balsam_job)

      return argojob

