#!/usr/bin/env python

import logging
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')
logger = logging.getLogger(__name__)
import optparse,os,sys,time,shutil,subprocess

import Computer
import JobTransaction
from AlpgenCmdFile import AlpgenCmdFile
from CondorJob import CondorJob,CondorJobStatus
import GridFtp
from update_unw_par_file import UpdateUnwParFile
sys.path.append('/users/hpcusers/svn/tools/python/athena')
from alpgen_unw_dataset import create_alpgen_dataset
from prun_gen_trf_grid import prun as prun_on_grid

GRIDFTP_PROTOCOL        = 'gsiftp://'
GRIDFTP_SERVERNAME      = 'atlasgridftp02.hep.anl.gov'
GRIDFTP_PATH            = '/grid/atlas/hpc/transfer'
ALPGEN_INSTALL          = '/users/hpcusers/svn/generators/alpgen/v214/x86/trunk/'
PYTHIA_EXE              = '/users/hpcusers/svn/generators/pythia/v8180/usercode/alpToHepmc/trunk/runPythiaOnAlpgen'
PYTHIA_INSTALL          = '/users/hpcusers/svn/generators/pythia/v8180/bgq/trunk'
PYTHIA_HEPMC_FILEBASE   = 'pythia_output'
PYTHIA_MERGE_HEPMC      = '/users/hpcusers/balsam/exe/merge_hepmc.py'
TO_HPC_PATH             = 'to_hpc'
FROM_HPC_PATH           = 'from_hpc'
PDF_FILENAME            = '/grid/atlas/hpc/common/alpgen/pdf/cteq6l1.tbl'
CONDOR_ENVIRONMENT      = 'LD_LIBRARY_PATH=/usr/lib64/openmpi/lib'

GRID_EVGEN_JOB_OPTS     = '/grid/atlas/hpc/common/athena/MC12JobOpts-00-11-59_v5_ANL.tar.gz'
GRID_SCRIPT             = '/grid/atlas/hpc/common/athena/gen_trf_grid.py'


ALPGEN_WARMUP_JOB_POLL_TIME_SEC     = 300
ALPGEN_WARMUP_JOB_MAX_WAIT_SEC      = 60*60*24
ALPGEN_WEIGHTED_JOB_POLL_TIME_SEC   = 300
ALPGEN_WEIGHTED_JOB_MAX_WAIT_SEC    = 60*60*24*7
ALPGEN_UNWEIGHT_JOB_POLL_TIME_SEC   = 300
ALPGEN_UNWEIGHT_JOB_MAX_WAIT_SEC    = 60*60*24
PYTHIA_JOB_POLL_TIME_SEC            = 300
PYTHIA_JOB_MAX_WAIT_SEC             = 60*60*24


def main():
   # print command
   logger.info( str(sys.argv) )
   # parse command line
   parser = optparse.OptionParser(description='job payload for running Pythia8',usage=sys.argv[0]+' -e <evts-per-iter> -i <num-iter> -n <num-weighted> -m <machine>\n         -p <process> -o <num-nodes> -c <cpus-per-node> -a <alpgen-input-file>\n        -t <walltime> [-s] [-u] [-b HPC_BINARY_NAME]\n   use [-h,--help] for details')
   parser.add_option('-e','--evts-per-iter',dest='evts_per_iter',help='number of events per warmup iteration',type='int')
   parser.add_option('-i','--num-iter',dest='numiters',help='number of iterations for the warmup',type='int')
   parser.add_option('-w','--warmup-weighted',dest='num_warmup',help='number of event to in the warmup, after the iterations complete',type='int')
   parser.add_option('-n','--num-weighted',dest='num_weighted',help='number of weighted events to generate.',type='int')
   parser.add_option('-m','--machine',dest='machine_name',help='specify the machine to which the job is submitted')
   parser.add_option('-p','--process',dest='process',help='define the process to generate, 2Q,4Q,hjet,top,wjet,zjet,Njet,etc.')
   parser.add_option('-o','--num-nodes',dest='numnodes',help='number of nodes to use on destination machine',type='int')
   parser.add_option('-c','--cpus-per-node',dest='cpus_per_node',help='number of CPUs per node to use on destination machine',type='int')
   parser.add_option('-a','--alpgen-input',dest='alpgen_input_file',help='The AlpGen input file which carries all the options for this generation job')
   parser.add_option('-u','--skip-hpc',dest='skip_hpc',help='Skip submission to HPC.',default=False,action='store_true')
   parser.add_option('-t','--wall-time',dest='walltime',help='The wall time to submit to the queue in minutes.',type='int')
   parser.add_option('-r','--resubmit',dest='resubmitjobid',help='Pass the job ID and this job will be resubmitted at the HPC weighted event gen step')
   parser.add_option('-b','--hpc-binary',dest='binary',help='name of binary to use on HPC (if using a custom)')
   parser.add_option('-j','--jobConfig',dest='jobConfig',help='Athena Generate_trf jobConfig setting')
   parser.add_option('-d','--inDS',dest='inDS',help='Dataset which is input to the Generate_trf call')
   parser.add_option('-f','--outDS',dest='outDS',help='Dataset which is output from the Generate_trf call')
   options,args = parser.parse_args()
   

   if options.numiters is None:
      parser.error('Must define the number of warmup iterations')
   if options.machine_name is None:
      parser.error('Must define the destination machine')
   if options.process is None:
      parser.error('Must define the process to generate')
   if not Computer.found_computer(options.machine_name):
      parser.error('Machine name '+options.machine_name+' not found.')
   if options.numnodes is None:
      parser.error('Must define the number of nodes to use')
   if options.cpus_per_node is None:
      parser.error('Must define the number of CPUs per node to use')
   if options.evts_per_iter is None:
      parser.error('Must define the number of events per warmup iteration')
   if options.num_weighted is None:
      parser.error('Must define the number of weighted events to produce')
   if options.num_warmup is None:
      parser.error('Must define the number of weighted events to produce in the warmup step.')
   if options.alpgen_input_file is None:
      parser.error('Must define the AlpGen input file')
   if options.walltime is None:
      parser.error('Must specify a wall time')
   if options.jobConfig is None:
      parser.error('Must specify -j')
   if options.inDS is None:
      parser.error('Must specify -d')
   if options.outDS is None:
      parser.error('Must specify -f')
   
   #options.numnodes = int(options.numnodes)
   #options.cpus_per_node = int(options.cpus_per_node)
   #options.numiters = int(options.numiters)
   #options.num_weighted = int(options.num_weighted)
   #options.evts_per_iter = int(options.evts_per_iter)

   # check that computer configuration is allowed
   if not Computer.check_computer_config(options.machine_name,options.numnodes,options.cpus_per_node):
      logger.error('ERROR in computer configuration.')
      return -1

   # derive alpgen exe from process name
   ALPGEN_EXE = ALPGEN_INSTALL + '/' + options.process + 'work/' + options.process + 'gen90'
   if not os.path.isfile(ALPGEN_EXE):
      logger.error('No alpgen executable found, path give:\n     ' + ALPGEN_EXE)
      return -2

   HPC_ALPGEN_EXE = os.path.basename(ALPGEN_EXE)
   if options.binary is not None:
      HPC_ALPGEN_EXE = options.binary

   # check pythia executable exists
   if not os.path.isfile(PYTHIA_EXE):
      logger.error('No pythia executable found, path given:\n    ' + PYTHIA_EXE)
      return -100
   
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
   if(user == 'apf'): # AutoPyFactory
      user= os.environ.get('prodUserID','nobody')
   jobID = taskID + '0'

   if options.resubmitjobid is not None:
      jobID = int(options.resubmitjobid)

   TOP_PATH = os.getcwd() # directory in which script was run
   RUNPATH = os.path.join(TOP_PATH,str(jobID)) # directory in which to store files
   if not os.path.exists(RUNPATH):
      os.makedirs(RUNPATH) # make directories recursively like 'mkdir -p'
   
   logger.info('JobID: ' + str(jobID))



   # move to run directory
   os.chdir(RUNPATH)

   # write an MD5 checksum for the local ALPGEN_EXE and PYTHIA_EXE
   write_checksum(ALPGEN_EXE,'checksum_alpgen.md5')
   write_checksum(PYTHIA_EXE,'checksum_pythia.md5')

   
   # add the top_dir to alpgen_input_filename, but take care that relative path could be included in name
   alpgen_input_filename = options.alpgen_input_file
   if len(alpgen_input_filename.split('/')) == 1: # string is just a filename, does not include a path
      # attach top_dir path
      alpgen_input_filename = os.path.join(TOP_PATH,alpgen_input_filename)
   elif len(alpgen_input_filename.split('/')) > 1: # string contains a path
      if alpgen_input_filename[0] == '.': # relative path was passed
         # add top dir and normalize to remove relative path
         alpgen_input_filename = os.path.normpath(os.path.join(TOP_PATH,alpgen_input_filename))

   # create gridFTP URLs 
   gridftp_from_hpc_path   = GRIDFTP_PROTOCOL + GRIDFTP_SERVERNAME + os.path.join(GRIDFTP_PATH,FROM_HPC_PATH,str(jobID)) + '/'
   gridftp_to_hpc_path     = GRIDFTP_PROTOCOL + GRIDFTP_SERVERNAME + os.path.join(GRIDFTP_PATH,TO_HPC_PATH,  str(jobID)) + '/'
  
    # create Alpgen input file
   input_file = AlpgenCmdFile()
   input_file.read(alpgen_input_filename)
   
   ################################################
   #####         run warmup grid generation
   ###############################################
   if options.resubmitjobid is None:
      alpgen_warmup_step(ALPGEN_EXE,
                         alpgen_input_filename,
                         options.evts_per_iter,
                         options.numiters,
                         options.num_warmup,
                        )
   
   #################################################
   #####        run weighted event generation
   ################################################
   
     
   # run weighted event generation
   if not options.skip_hpc:
      alpgen_weighted_event_gen(HPC_ALPGEN_EXE,
                                 alpgen_input_filename,
                                 options.num_weighted,
                                 options.numnodes,
                                 options.cpus_per_node,
                                 gridftp_to_hpc_path,
                                 jobID,
                                 user,
                                 options.machine_name,
                                 options.walltime,
                              )
   
      # copy files from GridFTP server back to current path
      GridFtp.globus_url_copy(gridftp_from_hpc_path + '/','./') 
      # fix event number in par file
      update_alpgen_weighted_par_file(alpgen_input_filename)
      

      ########################################################
      ########           run event unweighting 
      #######################################################
      alpgen_event_unweighting(ALPGEN_EXE,alpgen_input_filename)
   
      # add parameters to par file
      UpdateUnwParFile(input_file.filename_base + '_unw.par')



   ########################################################
   ########       send job off to grid for finishing
   #######################################################
   submit_generate_trf_on_grid(
                               alpgen_base=input_file.filename_base,
                               evts_per_file=5000,
                               inDS = options.inDS,
                               outDS = options.outDS,
                               ecmEnergy=8000,
                               runNumber=options.jobConfig.split('.')[1],
                               jobConfig = options.jobConfig,
                               evgenJobOpts = GRID_EVGEN_JOB_OPTS,
                               grid_script = GRID_SCRIPT,
                               )


def alpgen_warmup_step(alpgen_exe,alpgen_input_filename,events_per_iteration,num_iter,num_warmup):
   
   logger.info('Begin AlpGen warmup')

   # create Alpgen input file
   input_file = AlpgenCmdFile()
   input_file.read(alpgen_input_filename)
   input_file.imode        = 0
   input_file.start_with   = 0
   input_file.nevt         = events_per_iteration
   input_file.nitr         = num_iter
   input_file.last_nevt    = num_warmup
   
   custom_input_filename = 'input.0'
   input_file.write(custom_input_filename)
   # make link to input.0 named 'input' because alpgen expects a file called 'input'
   if os.path.exists('input'):
      os.remove('input')
   os.symlink(custom_input_filename,'input')

   condor_job = CondorJob(executable = alpgen_exe,
                          transfer_input_files  = PDF_FILENAME + ',input',
                          transfer_output_files = input_file.filename_base + '.grid1,' + input_file.filename_base + '.grid2',
                          output = 'warmup.stdout',
                          error  = 'warmup.stderr',
                          log    = 'warmup.condor_log',
                          environment = CONDOR_ENVIRONMENT,
                         )
   # write condor job file for record keeping
   condor_job.write("warmup.condor")
   condor_job.submit()
   
   start_time = time.time()
   while True:
      # see if job has completed
      status = condor_job.get_status()
      if status == CondorJobStatus.Completed:
         break
      elif status == CondorJobStatus.Holding:
         # job is held so it failed
         stdout,stderr = condor_job.analyze()
         if stdout is not None and stderr is not None:
            logger.error('Job is held. Here is the condor_q -analyze output:\nstdout = \n' + stdout + '\nstderr = \n' + stderr + '\n')
         condor_job.remove()
         logger.fatal('Failed to run Alpgen Warmup Step.')
         sys.exit(-1)
      
      # don't wait forever
      run_time = time.time() - start_time
      if run_time > ALPGEN_WARMUP_JOB_MAX_WAIT_SEC:
         logger.warning('CondorJob timed out')
         break
      # sleep and try again in a bit
      time.sleep(ALPGEN_WARMUP_JOB_POLL_TIME_SEC)

   # remove symlink
   os.remove('input')
   logger.info('Completed AlpGen Warmup')

def alpgen_weighted_event_gen(alpgen_exe,
                              alpgen_input_filename,
                              num_weighted_evts,
                              num_nodes,
                              cpus_per_node,
                              gridftp_to_hpc_path,
                              job_id,
                              user,
                              machine_name,
                              walltime
                             ):
   logger.info('Begin AlpGen weighted event generation')
   # parse input file
   input_file = AlpgenCmdFile()
   input_file.read(alpgen_input_filename)
   alpgen_basename = input_file.filename_base
   input_file.imode        = 1
   input_file.start_with   = 2
   input_file.nevt         = 0
   input_file.nitr         = 0
   input_file.last_nevt    = (num_weighted_evts/(num_nodes*cpus_per_node))
   custom_input_filename = 'input.1'
   input_file.write(custom_input_filename)
   # copy input file to gridftp area
   GridFtp.globus_url_copy(custom_input_filename,gridftp_to_hpc_path)
   # copy files Alpgen grid files, PDF to gridftp area for job submission
   GridFtp.globus_url_copy(alpgen_basename+'.grid1',gridftp_to_hpc_path)
   GridFtp.globus_url_copy(alpgen_basename+'.grid2',gridftp_to_hpc_path)
   GridFtp.globus_url_copy(PDF_FILENAME            ,gridftp_to_hpc_path)

   
   # build job description for weighted event generation
   jobDesc = JobTransaction.JobDescription(
                              job_id            = job_id,
                              user              = user,
                              num_evts          = num_weighted_evts,
                              exe               = os.path.basename(alpgen_exe),
                              nodes             = num_nodes,
                              cpus_per_node     = cpus_per_node,
                              task_duration_minutes = walltime,
                              input_files       = (alpgen_basename+'.grid1,'
                                                  +alpgen_basename+'.grid2,'
                                                  + PDF_FILENAME+','
                                                  +custom_input_filename
                                                  ),
                              output_files      = alpgen_basename+'.par,'+alpgen_basename+'.wgt',
                              postprocess       = 'postsubmit.sh',
                              preprocess        = 'presubmit.sh',
                              preprocess_args   = str(num_nodes*cpus_per_node) + ' ' + custom_input_filename,
                              transfer_out      = JobTransaction.TRUE,
                              transfer_in       = JobTransaction.TRUE
                           )
   #jobDesc.write_job_xml_to_file(PYTHIA_JOB_XML)
   
   msgq = JobTransaction.JobMsgHandler(job_id = jobDesc.job_id)
   msgq.send_submit_job_msg(
                            machine_name     = machine_name,
                            operation        = 'SUBMIT',
                            job_description  = jobDesc,
                           )
   
   start_time = time.time()
   while True:
      msg = msgq.recv_job_status_msg()
      if msg.contains(JobTransaction.BalsamMessage.SUCCEEDED):
         logger.debug('job succeeded')
         if msg.contains(JobTransaction.BalsamMessage.SUBMIT_DISABLED):
            logger.info('job succeeded, but HPC job submission is disabled so exiting')
            sys.exit(0)
         break
      elif msg.message is JobTransaction.BalsamMessage.NO_MESSAGE:
         logger.debug('waiting for message from machine: ' + machine_name)
         # don't wait forever
         run_time = time.time() - start_time
         if run_time > ALPGEN_WEIGHTED_JOB_MAX_WAIT_SEC:
            logger.warning('HPC Job timed out')
            break
         # sleep for a bit and try again
         time.sleep(ALPGEN_WEIGHTED_JOB_POLL_TIME_SEC)
      elif msg.contains(JobTransaction.BalsamMessage.FAILED):
         logger.error('job failed HPC submission, message = ' + str(msg))
         sys.exit(-1)

   logger.info('Completed AlpGen weighted event generation')



def alpgen_event_unweighting(alpgen_exe,alpgen_input_filename):
   logger.info('Begin AlpGen unweighting') 
   # create Alpgen input file
   input_file = AlpgenCmdFile()
   input_file.read(alpgen_input_filename)
   input_file.imode        = 2
   input_file.start_with   = 1
   input_file.nevt         = 0
   input_file.nitr         = 0
   input_file.last_nevt    = 999999999
   custom_input_filename = 'input.2'
   input_file.write(custom_input_filename)
   # create symlink for alpgen
   if os.path.exists('input'):
      os.remove('input')
   os.symlink(custom_input_filename,'input')
      
   condor_job = CondorJob(executable = alpgen_exe,
                          transfer_output_files = (input_file.filename_base + '.unw,'
                                                  +input_file.filename_base + '_unw.par'
                                                  ),
                          transfer_input_files  = (input_file.filename_base + '.wgt,'
                                                  +input_file.filename_base + '.par,'
                                                  +PDF_FILENAME + ','
                                                  +'input'
                                                  ),
                          output       = 'unweight.stdout',
                          error        = 'unweight.stderr',
                          log          = 'unweight.condor_log',
                         )
   # write condor job file for record keeping
   condor_job.write("unweight.condor")
   condor_job.submit()

   start_time = time.time()
   while True:
      # see if job has completed
      status = condor_job.get_status()
      if status == CondorJobStatus.Completed:
         break
      elif status == CondorJobStatus.Holding:
         # job is held so it failed
         stdout,stderr = condor_job.analyze()
         if stdout is not None and stderr is not None:
            logger.error('Job is held. Here is the condor_q -analyze output:\nstdout = \n' + stdout + '\nstderr = \n' + stderr + '\n')
         condor_job.remove()
         logger.fatal('Failed to run Alpgen Unweighting Step.')
         sys.exit(-1)

      
      # don't wait forever
      run_time = time.time() - start_time
      if run_time > ALPGEN_UNWEIGHT_JOB_MAX_WAIT_SEC:
         logger.warning('CondorJob timed out')
         break
      # sleep and try again in a bit
      time.sleep(ALPGEN_UNWEIGHT_JOB_POLL_TIME_SEC)

   # remove symlink
   os.remove('input')
   logger.info('Completed AlpGen unweighting')
                          


def shower_alpgen_with_pythia(pythia_exe,alpgen_input_filename):

   logger.info('Begin Pythia showering of events')

   input_file = AlpgenCmdFile()
   input_file.read(alpgen_input_filename)
   
   # build arguments for pythia executable
   args  = ' -x ' + os.path.join(PYTHIA_INSTALL,'xmldoc')
   args += ' -o ' + PYTHIA_HEPMC_FILEBASE + '.hepmc'
   args += ' -a ' + input_file.filename_base


   condor_job = CondorJob(executable = PYTHIA_EXE,
                          transfer_input_files     = (input_file.filename_base + '.unw,'
                                                     +input_file.filename_base + '_unw.par,'
                                                     +PYTHIA_MERGE_HEPMC
                                                     ),
                          transfer_output_files    = PYTHIA_HEPMC_FILEBASE + '.hepmc',
                          postcmd                  = 'merge_hepmc.py',
                          output                   = 'pythia.stdout',
                          error                    = 'pythia.stderr',
                          log                      = 'pythia.condor_log',
                          arguments                = args,
                         )
   # write condor job file for record keeping
   condor_job.write("pythia.condor")
   condor_job.submit()

   start_time = time.time()
   while True:
      # see if job has completed
      status = condor_job.get_status()
      if status == CondorJobStatus.Completed:
         break
      elif status == CondorJobStatus.Holding:
         # job is held so it failed
         stdout,stderr = condor_job.analyze()
         if stdout is not None and stderr is not None:
            logger.error('Job is held. Here is the condor_q -analyze output:\nstdout = \n' + stdout + '\nstderr = \n' + stderr + '\n')
         condor_job.remove()
         logger.fatal('Failed to run Pythia showering step.')
         sys.exit(-1)

      
      # don't wait forever
      run_time = time.time() - start_time
      if run_time > PYTHIA_JOB_MAX_WAIT_SEC:
         logger.warning('CondorJob timed out')
         break
      # sleep and try again in a bit
      time.sleep(PYTHIA_JOB_POLL_TIME_SEC)
   logger.info('Completed Pythia showering')

def submit_generate_trf_on_grid(
                               alpgen_base,
                               evts_per_file,
                               inDS,
                               outDS,
                               ecmEnergy,
                               runNumber,
                               jobConfig,
                               evgenJobOpts,
                               grid_script,
                               site='MWT2_UC_SCRATCHDISK',
                               asetup='17.2.11.12,AtlasProduction',
                               randomSeed=12345678,
                               firstEvent=1,
                               outputEVNTFile='myOutputEVNTFile.root',
                               ):
   # first divide up the alpgen file
   try:
      create_alpgen_dataset(alpgen_base,evts_per_file,inDS,site,asetup)
   except:
      logger.error(' exception caught while creating and uploading alpgen dataset: ' + str(sys.exc_info()[1]))
      raise
   # second submit job with prun
   try:
      prun_on_grid(
           asetup,
           ecmEnergy,
           randomSeed,
           outputEVNTFile,
           runNumber,
           firstEvent,
           jobConfig,
           evgenJobOpts,
           inDS,
           outDS,
           grid_script,
          )
   except:
      logger.error(' exception caught while submitting prun job: ' + str(sys.exc_info()[1]))
      raise

           
           
   

def update_alpgen_weighted_par_file(alpgen_input_filename):
   
   input_file = AlpgenCmdFile()
   input_file.read(alpgen_input_filename)
   filename_base = input_file.filename_base
   
   par_filename = filename_base + '.par'
   weighted_events_filename = filename_base + '.wgt'
   new_par_filename = filename_base + '.par.new'

   # read the number of events in the weighted event file
   p = subprocess.Popen(['wc','-l',weighted_events_filename],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   stdout,stderr = p.communicate()
   if p.returncode != 0:
      logger.error('word count command failed. stdout: \n' + stdout + '\n stderr: \n' + stderr + '\n')
   new_wgt_event_count = float(stdout.split(' ')[0])


   par_file = open(par_filename,'r')
   new_par_file = open(new_par_filename,'w')
   
   # find the line where the number of weighted events are located
   numbers_are_next = False
   for line in par_file:
      if line.find('number wgted evts in the file') >= 0:
         # write current line to new file
         numbers_are_next = True
      elif numbers_are_next:
         # get next line with numbers on it, format is "# wgt events"  "sigma" "error" "maxwgt"
         number_strings = line.split()
         print len(number_strings), str(number_strings)
         sigma = number_strings[1]
         if sigma[0] is '.':
            sigma = '0' + sigma
         error = number_strings[2]
         if error[0] is '.':
            error = '0' + error
         maxwgt = number_strings[3]
         if maxwgt[0] is '.':
            maxwgt = '0' + maxwgt
         new_line = ("%15.1f %s %s %s\n" % (new_wgt_event_count,sigma,error,maxwgt))
         logger.debug('old par line: ' + line)
         logger.debug('new par line: ' + new_line)
         new_par_file.write(new_line)
         numbers_are_next = False
         continue
      new_par_file.write(line)

   par_file.close()
   new_par_file.close()

   os.remove(par_filename)
   shutil.copyfile(new_par_filename,par_filename)
         
def write_checksum(exe,filename = None):
   p = subprocess.Popen(['openssl','md5',exe],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   stdout,stderr = p.communicate()
   p = None

   if filename is not None:
      file = open(filename,'w')
      file.write(stdout)
      file.close()
      file = None
   return stdout



if __name__ == '__main__':
   #update_alpgen_weighted_par_file('alpout')
   main()



