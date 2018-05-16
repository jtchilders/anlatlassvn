#!/usr/bin/env python
import classad,htcondor
import logging,subprocess,optparse

logger = logging.getLogger(__name__)

def main():
   parser = optparse.OptionParser(description='Submit a job to condor via the command line or output a condor job file.')
   parser.add_option('-u','--universe',dest='universe',help='Set the Condor universe.',default='vanilla')
   parser.add_option('-e','--executable',dest='executable',help='The executable.')
   parser.add_option('-a','--arguments',dest='arguments',help='Set the executable arguments.')
   parser.add_option('-r','--requirements',dest='requirements',help='Set the requirements on the target machines.')
   parser.add_option('-g','--getenv',dest='getenv',help='Set the Condor environment from the current shell. Can only be used if writing the condor job out to a submit file because job submission in this script is done via a sub-shell.',action='store_true',default=False)
   parser.add_option('-n','--environment',dest='environment',help='Environment variables to set in job shell. Use quotations here and follow the new condor format.')
   parser.add_option('-i','--input',dest='input',help='File to be used to satisfy stdin if needed.')
   parser.add_option('-o','--output',dest='output',help='File for stdout log.',default='stdout.txt')
   parser.add_option('-b','--error',dest='error',help='File for stderr log.',default='stderr.txt')
   parser.add_option('-l','--log',dest='condor_log',help='File for condor log.',default='condor_log.txt')
   parser.add_option('-t','--transfer-off',dest='transfer_off',help='Turns off file transfers for jobs',action='store_true',default=False)
   parser.add_option('-x','--transfer-when',dest='transfer_when',help='When should output files be transfer? ON_EXIT or ON_EXIT_OR_EVICT',default='ON_EXIT_OR_EVICT')
   parser.add_option('-c','--transfer-in',dest='transfer_in',help='Files to be transferred into job. Comma separated list.')
   parser.add_option('-d','--transfer-out',dest='transfer_out',help='Files to be transferred out of job. Comma separated list.')
   parser.add_option('-f','--precmd',dest='precmd',help='Command to execute before primary job.')
   parser.add_option('-j','--preargs',dest='preargs',help='Arguments to be passed to command that is executed before primary job.')
   parser.add_option('-k','--postcmd',dest='postcmd',help='Command to execute after primary job.')
   parser.add_option('-m','--postargs',dest='postargs',help='Arguments to be passed to command that is executed after primary job.')
   parser.add_option('-q','--queue',dest='queue',help='Queue name to which the job will be submitted.')
   parser.add_option('-w','--write-condor-file',dest='write_condor_file',help='If used, the Condor Job submit file will be written to the file specified.')
   parser.add_option('-s','--submit',dest='submit',help='Submit the job directly.',action='store_true',default=False)
   
   options,args = parser.parse_args()

   if options.executable is None:
      parser.error('Must specify -e')

   transfer_files = 'YES'
   if options.transfer_off:
      transfer_files = 'NO'

   
   job = CondorJob(
                universe                  = options.universe,
                executable                = options.executable,
                arguments                 = options.arguments,
                requirements              = options.requirements,
                getenv                    = options.getenv,
                environment               = options.environment,
                input                     = options.input,
                output                    = options.output,
                error                     = options.error,
                log                       = options.condor_log,
                should_transfer_files     = transfer_files,
                when_to_transfer_output   = options.transfer_when,
                transfer_input_files      = options.transfer_in,
                transfer_output_files     = options.transfer_out,
                precmd                    = options.precmd,
                preargs                   = options.preargs,
                postcmd                   = options.postcmd,
                postargs                  = options.postargs,
                queue                     = options.queue,
              )
   
   if options.write_condor_file is not None:
      job.write(options.write_condor_file)

   if options.submit:
      job.submit()

class CondorJobStatus:
   """ Status Codes for Condor Jobs """
   Idle        = 1
   Running     = 2
   Removed     = 3
   Completed   = 4
   Holding     = 5

   Names = [
            'NoStatusCode',
            'Idle',
            'Running',
            'Removed',
            'Completed',
            'Holding',
           ]

   Map2JobState = {
         # hard-code these for now to avoid circular import
                   #Idle:       models.Job.CREATED,
                   #Running:    models.Job.RUNNING,
                   #Removed:    models.Job.FAILED,
                   #Completed:  models.Job.FINISHED,
                   #Holding:    models.Job.HOLDING,
                   Idle:       'C',
                   Running:    'R',
                   Removed:    'F',
                   Completed:  'JF',
                   Holding:    'HD'
                  }

class CondorJob:
   ''' create a cmd file for a condor job '''

   def __init__(self,
                universe                  = 'vanilla',
                executable                = 'echo',
                arguments                 = None,
                requirements              = None,
                getenv                    = False,
                environment               = None,
                input                     = None,
                output                    = 'out.txt',
                error                     = 'err.txt',
                log                       = 'log.txt',
                should_transfer_files     = 'YES',
                when_to_transfer_output   = 'ON_EXIT',
                transfer_input_files      = None,
                transfer_output_files     = None,
                transfer_output_remaps    = None,
                precmd                    = None,
                preargs                   = None,
                postcmd                   = None,
                postargs                  = None,
                queue                     = None,
                stream_output             = True,
                stream_error              = True,
               ):
      self.universe                 = universe
      self.executable               = executable
      self.arguments                = arguments
      self.requirements             = requirements
      self.getenv                   = getenv
      self.environment              = environment
      self.input                    = input
      self.output                   = output
      self.error                    = error
      self.log                      = log
      self.should_transfer_files    = should_transfer_files
      self.when_to_transfer_output  = when_to_transfer_output
      self.transfer_input_files     = transfer_input_files
      self.transfer_output_files    = transfer_output_files
      self.transfer_output_remaps   = transfer_output_remaps
      self.precmd                   = precmd
      self.preargs                  = preargs
      self.postcmd                  = postcmd
      self.postargs                 = postargs
      self.queue                    = queue
      self.stream_output            = stream_output
      self.stream_error             = stream_error

      self.ad           = None
      self.submitted_ad = None
      self.clusterId    = None
      
   def __str__(self):
      
      txt = ''
      if self.universe:                 txt  = 'Universe                 = ' + self.universe + '\n'
      if self.arguments:                txt += 'Arguments                = ' + self.arguments + '\n'
      if self.requirements:             txt += 'Requirements             = ' + self.requirements + '\n'
      if self.getenv:                   txt += 'GetEnv                   = ' + self.getenv + '\n'
      if self.input:                    txt += 'Input                    = ' + self.input + '\n'
      if self.output:                   txt += 'Output                   = ' + self.output + '\n'
      if self.error:                    txt += 'Error                    = ' + self.error + '\n'
      if self.log:                      txt += 'Log                      = ' + self.log + '\n'
      if self.should_transfer_files:    txt += 'should_transfer_files    = ' + self.should_transfer_files + '\n'
      if self.when_to_transfer_output:  txt += 'when_to_transfer_output  = ' + self.when_to_transfer_output + '\n'
      if self.executable:               txt += 'Executable               = ' + self.executable + '\n'
      if self.precmd:                   txt += '+PreCmd                  = "' + self.precmd + '"\n'
      if self.preargs:                  txt += '+PreArguments            = "' + self.preargs + '"\n'
      if self.postcmd:                  txt += '+PostCmd                 = "' + self.postcmd + '"\n'
      if self.postargs:                 txt += '+PostArguments           = "' + self.postargs + '"\n'
      if self.transfer_input_files:     txt += 'transfer_input_files     = ' + self.transfer_input_files + '\n'
      if self.transfer_output_files:    txt += 'transfer_output_files    = ' + self.transfer_output_files + '\n'
      if self.transfer_output_remaps:   txt += 'transfer_output_remaps   = "' + ';'.join(str(a)+' '+str(b) for a,b in self.transfer_output_remaps.iteritems()) + '"\n'
      if self.environment:              txt += 'environment              = ' + self.environment + '\n'
      txt                                   += 'StreamOut                = ' + str(self.stream_output) + '\n'
      txt                                   += 'StreamErr                = ' + str(self.stream_error) + '\n'
      if self.queue:                    txt += 'Queue ' + self.queue + '\n'
      else:                             txt += 'Queue\n'

      return txt
   
   def write(self,filename):

      try:
         file = open(filename,'w')
      except IOError,e:
         logger.error('CondorJob.write: Failed to open file ' + filename + ', exception: ' + str(e) )
         print -1

      try:
         file.write(str(self))
      except IOError,e:
         logger.error('CondorJob.write: Failed to write to file ' + filename + ', exception: ' + str(e) )
         return -2

      file.close()
      return 0

   def create_classAd(self):

      # create empty class ad
      self.ad = classad.ClassAd()
      # and fill it with the parameters defined for this job
      if self.universe:                 self.ad['Universe']             = self.universe
      if self.executable:               self.ad['Cmd']                  = self.executable
      if self.arguments:                self.ad['Arguments']            = self.arguments
      if self.requirements:             self.ad['Requirements']         = self.requirements
      if self.environment:              self.ad['Env']                  = self.environment
      if self.input:                    self.ad['In']                   = self.input
      if self.output:                   self.ad['Out']                  = self.output
      if self.error:                    self.ad['Err']                  = self.error
      if self.log:                      self.ad['UserLog']              = self.log
      if self.should_transfer_files:    self.ad['ShouldTransferFiles']  = self.should_transfer_files
      if self.when_to_transfer_output:  self.ad['WhenToTransferOutput'] = self.when_to_transfer_output
      if self.transfer_input_files:     self.ad['TransferInput']        = self.transfer_input_files
      if self.transfer_output_files:    self.ad['TransferOutput']       = self.transfer_output_files
      if self.transfer_output_remaps:   self.ad['TransferOutputRemaps'] = ';'.join(str(a)+' '+str(b) for a,b in self.transfer_output_remaps.iteritems())
      if self.precmd:                   self.ad['PreCmd']               = self.precmd
      if self.preargs:                  self.ad['PreArguments']         = self.preargs
      if self.postcmd:                  self.ad['PostCmd']              = self.postcmd
      if self.postargs:                 self.ad['PostArguments']        = self.postargs
      self.ad['stream_error']                                           = self.stream_error
      self.ad['stream_output']                                          = self.stream_output

      logger.debug(str(self.ad))
      return self.ad


   def submit(self):
      
      if self.ad is None:
         self.create_classAd()
      
      # get the scheduler
      schedd = htcondor.Schedd()
      # empty list will be passed to submit function and filled with an ad for the submited job
      submitted_ads = []
      cluster_id = schedd.submit(self.ad,1,False,submitted_ads)
      if cluster_id > 0:
         # job submitted correctly
         logger.debug('CondorJob.submit_ad: Job Submitted, Cluster ID: ' + str(cluster_id))
         # only submitted one job so list should have one entry
         self.submitted_ad = submitted_ads[0]
         self.clusterId = self.submitted_ad['ClusterId']
      else:
         logger.debug('CondorJob.submit_ad: Job Submission Failed, cluster_id = ' + str(cluster_id))
         return -1

      return 0

   def get_status(self):
      # some consistency checks
      if self.submitted_ad is None:
         logger.error('CondorJob.get_status: No job submitted to monitor')
         return -1
      if self.clusterId is None:
         logger.error('CondorJob.get_status: No job cluster id, job likely not submitted yet')
         return -2

      # get scheduler
      schedd = htcondor.Schedd()
      
      # get list of jobs with this jobs cluster ID (should only be one)
      try:
         status_jobs = schedd.query('ClusterId == ' + str(self.clusterId))
      except IOError,e:
         logger.debug('CondorJob.get_status: No job returned by query call, usually means the job is finished. Exception: ' + str(e))
      
      # should only be one job
      manyjobs = False
      if len(status_jobs) == 0:
         # if list is empty the job has completed or killed
         # so check in the history for the last state
         if(self.job_finished_successfully()):
            return CondorJobStatus.Completed
         else:
            logger.debug('CondorJob.get_status: ad list is empty, but check of history failed.')
            return CondorJobStatus.Removed
      if len(status_jobs) > 1:
         # retrieved more than one job (should not reach this)
         logger.debug('CondorScheduler.get_status: Retrieved more than one statuses, should not possible since the clusterID and processID uniquely define a condor job. Will only treat the first job.')
         manyjobs = True

      # use only the first one in the list
      status_ad = status_jobs[0]
      
      state = status_ad['JobStatus']

      logger.debug('CondorJob.get_status: job status = ' + CondorJobStatus.Names[state])

      if state == CondorJobStatus.Holding:
         logger.debug('CondorJob.get_status: stdout = \n' )

      return state

   def job_finished_successfully(self):
      condor_history_exe = 'condor_history'
      if self.clusterId:
         clusterId = str(self.clusterId)

         # call condor_history
         try:
            p = subprocess.Popen([condor_history_exe,clusterId],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
         except OSError,e:
            logger.exception('CondorJob.get_job_finished_state: Error calling condor_history'+ str(e))
         except ValueError,e:
            logger.exception('CondorJob.get_job_finished_state: Error calling condor_history'+ str(e))

         # get stdout and stderr
         out,err = p.communicate()
         lines = out.split('\n')
         if len(lines) > 3:
            logger.info('CondorJob.get_job_finished_state: Something strange is going on. Retrieved more than one history:\n'+out)
         # only care about second line and the first word should be the cluster id
         line = lines[1]
         words = line.split()
         history_cluster_id = words[0].split('.')[0]
         if int(history_cluster_id) != int(clusterId):
            logger.info('CondorJob.get_job_finished_state: Something is strange. The cluster id is mismatched, from job: '+clusterID+'; from condor_history: '+history_cluster_id)
            return False
         history_job_finished_state = words[5]

         if history_job_finished_state == 'C':
            return True

         return False

   def remove(self):
      if self.clusterId is not None:
         condor_rm_exe = 'condor_rm'

         p = subprocess.Popen([condor_rm_exe,str(self.clusterId)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
         out,err = p.communicate()

         if p.returncode is not 0:
            logger.warning('CondorJob.remove: Failed to remove job with cluster id = ' + str(self.clusterId) 
                           + '.\n  stdout = \n' + out + '\n stderr = \n' + err + '\n'
                          )
   
   def analyze(self):
      if self.clusterId is not None:
         condor_q_exe = 'condor_q'

         p = subprocess.Popen([condor_q_exe,'-analyze',str(self.clusterId)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
         out,err = p.communicate()

         if p.returncode is not 0:
            logger.warning('CondorJob.analyze: Failed to analyze job with cluster id = ' + str(self.clusterId) 
                           + '.\n  stdout = \n' + out + '\n stderr = \n' + err + '\n'
                          )

         return out,err
      return None,None





if __name__ == '__main__':
   logging.basicConfig(level=logging.INFO)
   main()
