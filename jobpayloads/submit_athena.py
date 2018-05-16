#!/usr/bin/env python
import os,sys,optparse,logging,time
import ArgoJobs.GridFtp as GridFtp
sys.path.append('/users/hpcusers/balsam/balsam_deploy/common_core')
from MessageInterface import MessageInterface
from ArgoJob import ArgoJob
from BalsamJob import BalsamJob
logger = logging.getLogger(__name__)

GRIDFTP_SERVER = 'gsiftp://atlasgridftp02.hep.anl.gov'
JOB_PATH = '/grid/atlas/hpc/argo/job'
ATHENA_CONDOR_SCRIPT = 'athena_condor_job.py'

def main():
   logging.basicConfig(level=logging.INFO)

   parser = optparse.OptionParser(description='')
   parser.add_option('-g','--run-gen-trf',dest='run_gen_trf',help='Run Generate_trf instead of athena.',action='store_true',default=False)
   parser.add_option('-a','--run-atlasG4-trf',dest='run_atlasG4_trf',help='Run AtlasG4_trf instead of athena.',action='store_true',default=False)
   parser.add_option('-d','--run-digi-trf',dest='run_digi_trf',help='Run Digi_trf instead of athena.',action='store_true',default=False)
   parser.add_option('-r','--run-reco-trf',dest='run_reco_trf',help='Run Reco_trf instead of athena.',action='store_true',default=False)
   parser.add_option('-y','--args',dest='args',help='Arguments to pass to athena or trf.',default='')
   parser.add_option('-j','--asetup-args',dest='asetup_args',help='Arguments used to setup Athena environment',default='')
   parser.add_option('-k','--atlas-local',dest='atlas_local_root_base',help='Path to ATLAS Local setup, typically called ATLAS_LOCAL_ROOT_BASE.',default='/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase')
   parser.add_option('-i','--input-files',dest='input_files',help='Comma separated list of files to transfer into Condor job',default='')
   parser.add_option('-o','--output-files',dest='output_files',help='Comma separated list of files to transfer out of Conodr job',default='')
   parser.add_option('-s','--site',dest='site',help='Target Balsam Site',default='argo_cluster')
   parser.add_option('-n','--num-nodes',dest='numnodes',help='number of nodes to use on destination machine',type='int')
   parser.add_option('-c','--cpus-per-node',dest='cpus_per_node',help='number of CPUs per node to use on destination machine',type='int')
   parser.add_option('-t','--wall-time',dest='walltime',help='The wall time to submit to the queue in minutes.',type='int',default=60)
   parser.add_option('-x','--no-submit',dest='submit',help='do not submit the message to ARGO. For testing purposes.',action='store_false',default=True)
   options,args = parser.parse_args()
   
   if (options.numnodes is None or 
       options.cpus_per_node is None):
      parser.error('Options -n and -c are required.')

   trf_count = 0
   executable = ''
   if options.run_gen_trf:
      executable = '-g'
      trf_count += 1
   if options.run_atlasG4_trf:
      executable = '-a'
      trf_count += 1
   if options.run_digi_trf:
      executable = '-d'
      trf_count += 1
   if options.run_reco_trf:
      executable = '-r'
      trf_count += 1
   if trf_count > 1:
      parser.error(' options -g -a -d -r are mutually exclusive, only one can be specified. ')
   
   submit_athena(executable,
                  options.args,
                  options.numnodes,
                  options.cpus_per_node,
                  options.asetup_args,
                  options.atlas_local_root_base,
                  options.input_files,
                  options.output_files,
                  options.site,
                  options.submit,
                  options.walltime,
                 )

def submit_athena(executable,
                  args,
                  numnodes,
                  cpus_per_node,
                  asetup_args,
                  atlas_local_root_base = '/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase',
                  input_files = '',
                  output_files = '',
                  site = 'argo_cluster',
                  submit = True,
                  walltime = 60,
                 ):

   jobID=str(int(time.time()*1000000))
   logger.info('JobID: ' + str(jobID))
   user = os.environ.get('USER','nobody')
   if(user == 'apf'): # AutoPyFactory
      user= os.environ.get('prodUserID','nobody')

   input_url = GRIDFTP_SERVER + JOB_PATH + '/' + jobID
   output_url = input_url

   # copy input files to input url
   if submit:
      for file in input_files.split(','):
         GridFtp.globus_url_copy(file,input_url+'/')
   input_files = [os.path.basename(file) for file in input_files.split(',')]
   input_files_str = (str(input_files)).replace(' ','').replace('[','').replace(']','').replace("'",'').replace('"','')

   # build argument string for athena condor script
   exe_args  = executable
   exe_args += ' -k ' + atlas_local_root_base
   if len(args)         > 0: exe_args += ' -x "' + args + '" '
   if len(asetup_args)  > 0: exe_args += ' -j "' + asetup_args + '" '
   if len(input_files)  > 0: exe_args += ' -i '  + input_files_str
   if len(output_files) > 0: exe_args += ' -o '  + output_files

   logger.info(' argument string: ' + exe_args )

   # create ArgoJob
   job = ArgoJob()
   # create Balsam Subjob
   subjob = BalsamJob(
                      executable          = ATHENA_CONDOR_SCRIPT,
                      executable_args     = exe_args,
                      input_files         = input_files,
                      output_files        = output_files.split(','),
                      nodes               = numnodes,
                      processes_per_node  = cpus_per_node,
                      target_site         = site,
                      wall_minutes        = walltime,
                     )
   
   job.input_url = input_url
   job.output_url = output_url
   job.username = user
   job.email = 'jchilders@anl.gov'
   job.add_job(subjob)

   logger.info(' job text: \n' + job.serialize())

   # submit job if requested
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
      mi.send_msg(job.serialize(),routing_key)
      print ' done sending'
      mi.close()
      print ' closing connection'
   else:
      logger.info(' not submitting job ')



   
   



if __name__ == "__main__":
   main()
