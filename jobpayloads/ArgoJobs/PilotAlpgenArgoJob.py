#!/usr/bin/env python
import os,sys,logging,shutil,time,optparse,subprocess,imp
import urllib
logger = logging.getLogger(__name__)
from AlpgenInputFile import AlpgenInputFile
import GridFtp
sys.path.append('/users/hpcusers/argobalsam/production/argobalsam/common_core')
from ArgoJob import ArgoJob
from BalsamJob import BalsamJob

# Path where job specific directories can be created
JOB_DIRECTORY        = '/grid/atlas/hpc/argo/jobs'
GRID_FTP_PROTOCOL    = ''
EVGEN_JOB_OPTS_URL   = 'http://atlas-computing.web.cern.ch/atlas-computing/links/kitsDirectory/EvgenJobOpts/'
EVGEN_JOB_OPTS_URL   = 'http://childers.web.cern.ch/childers/'
PDF_FILENAME         = '/grid/atlas/hpc/common/alpgen/pdf/cteq6l1.tbl'

gridgen_rates = {
   "wcjet + 0 jets": {
      "mean": 42817.067769100002,
      "n": 1,
      "sigma": 0.0
   },
   "wcjet + 1 jets": {
      "mean": 28708.133971300002,
      "n": 1,
      "sigma": 0.0
   },
   "wcjet + 2 jets": {
      "mean": 16040.3299725,
      "n": 1,
      "sigma": 0.0
   },
   "wcjet + 3 jets": {
      "mean": 10243.569314799999,
      "n": 1,
      "sigma": 0.0
   },
   "wcjet + 4 jets": {
      "mean": 5285.8261550500001,
      "n": 1,
      "sigma": 0.0
   },
   "wjet + 0 jets": {
      "mean": 60841.301736166664,
      "n": 3,
      "sigma": 7502.6467456600003
   },
   "wjet + 1 jets": {
      "mean": 32770.894687366665,
      "n": 3,
      "sigma": 3192.9826627799998
   },
   "wjet + 2 jets": {
      "mean": 22507.118883725001,
      "n": 4,
      "sigma": 4007.9773398799998
   },
   "wjet + 3 jets": {
      "mean": 13561.274375000001,
      "n": 4,
      "sigma": 2390.4127591699998
   },
   "wjet + 4 jets": {
      "mean": 7210.8120002733331,
      "n": 3,
      "sigma": 241.59335965099999
   },
   "wjet + 5 jets": {
      "mean": 3703.1342722666668,
      "n": 3,
      "sigma": 151.680436114
   },
   "wjet + 6 jets": {
      "mean": 2991.5041282799998,
      "n": 1,
      "sigma": 0.0
   },
   "wqq + 0 jets": {
      "mean": 33590.217570499997,
      "n": 2,
      "sigma": 308.08751419999999
   },
   "wqq + 1 jets": {
      "mean": 19792.173689349998,
      "n": 2,
      "sigma": 439.96346325000002
   },
   "wqq + 2 jets": {
      "mean": 13229.046291049999,
      "n": 2,
      "sigma": 178.47532855
   },
   "wqq + 3 jets": {
      "mean": 7373.8153515650001,
      "n": 2,
      "sigma": 93.563991305000002
   },
   "zjet + 0 jets": {
      "mean": 63112.920989254548,
      "n": 11,
      "sigma": 17050.574052399999
   },
   "zjet + 1 jets": {
      "mean": 26265.810577933335,
      "n": 3,
      "sigma": 5654.4957092000004
   },
   "zjet + 2 jets": {
      "mean": 19864.982661524999,
      "n": 8,
      "sigma": 2577.9874545299999
   },
   "zjet + 3 jets": {
      "mean": 12080.364631363333,
      "n": 3,
      "sigma": 2443.5158631899999
   },
   "zjet + 4 jets": {
      "mean": 5172.4610842833335,
      "n": 3,
      "sigma": 1044.7715383100001
   },
   "zjet + 5 jets": {
      "mean": 3401.50323501,
      "n": 3,
      "sigma": 722.77086425699997
   }
}

evgen_rates = {
   "wcjet + 1 jets": {
      "mean": 1170.5685618699999,
      "n": 1,
      "sigma": 0.0
   },
   "wcjet + 2 jets": {
      "mean": 1026.39296188,
      "n": 1,
      "sigma": 0.0
   },
   "wcjet + 3 jets": {
      "mean": 752.68817204300001,
      "n": 1,
      "sigma": 0.0
   },
   "wcjet + 4 jets": {
      "mean": 427.35042735000002,
      "n": 1,
      "sigma": 0.0
   },
   "wjet + 0 jets": {
      "mean": 504.34435234199998,
      "n": 3,
      "sigma": 24.228758177900001
   },
   "wjet + 1 jets": {
      "mean": 727.42918124955554,
      "n": 9,
      "sigma": 270.15562179400001
   },
   "wjet + 2 jets": {
      "mean": 887.76875175787495,
      "n": 16,
      "sigma": 159.088415895
   },
   "wjet + 3 jets": {
      "mean": 597.69705879935998,
      "n": 50,
      "sigma": 124.89395675999999
   },
   "wjet + 4 jets": {
      "mean": 531.97867806187628,
      "n": 59,
      "sigma": 121.472154512
   },
   "wjet + 5 jets": {
      "mean": 399.2158004500107,
      "n": 280,
      "sigma": 20.474836054899999
   },
   "wjet + 6 jets": {
      "mean": 222.024866785,
      "n": 1,
      "sigma": 0.0
   },
   "wqq + 0 jets": {
      "mean": 1037.4161573375,
      "n": 2,
      "sigma": 59.7625260525
   },
   "wqq + 1 jets": {
      "mean": 1140.5007801449999,
      "n": 2,
      "sigma": 22.289917525
   },
   "wqq + 2 jets": {
      "mean": 885.35192371142853,
      "n": 7,
      "sigma": 14.5283281114
   },
   "wqq + 3 jets": {
      "mean": 555.00301398644115,
      "n": 34,
      "sigma": 6.14490852613
   },
   "zjet + 0 jets": {
      "mean": 350.68371436585716,
      "n": 7,
      "sigma": 227.15929025599999
   },
   "zjet + 1 jets": {
      "mean": 253.70941517899999,
      "n": 3,
      "sigma": 3.3838742758799998
   },
   "zjet + 2 jets": {
      "mean": 635.15884400699997,
      "n": 12,
      "sigma": 323.23087095800003
   },
   "zjet + 3 jets": {
      "mean": 108.45165891466667,
      "n": 3,
      "sigma": 1.1100724078099999
   },
   "zjet + 4 jets": {
      "mean": 371.3474725285846,
      "n": 13,
      "sigma": 132.45901845
   },
   "zjet + 5 jets": {
      "mean": 284.74205536867657,
      "n": 47,
      "sigma": 47.028652389900003
   }
}

pythia_shower_eff = {
   0: 0.74,
   1: 0.30,
   2: 0.16,
   3: 0.07,
   4: 0.03,
   5: 0.02,
}

class PilotAlpgenArgoJob:

   def __init__(self,
                ecm,
                run_number,
                job_config,
                evgen_job_opts,
                output_evnt_file,
                ami_tag,
                steering,
                nevts,
                username      = os.environ['USER'],
                gridgen_site  = 'argo_cluster',
                evgen_site    = 'argo_cluster',
                shower_site   = 'argo_cluster',
                job_directory = JOB_DIRECTORY,
               ):
      self.ecm                = ecm
      self.run_number         = run_number
      self.job_config         = job_config
      self.evgen_job_opts     = evgen_job_opts
      self.output_evnt_file   = output_evnt_file
      self.ami_tag            = ami_tag
      self.steering           = steering
      self.nevts              = nevts
      self.username           = username
      self.gridgen_site       = gridgen_site
      self.evgen_site         = evgen_site
      self.shower_site        = shower_site
      self.job_directory      = job_directory



   def get_argo_job(self):

      job_path                = get_unique_job_path(self.job_directory)
      logger.debug('job path: ' + job_path)

      ##----------------------------------------------------
      #   Create Alpgen Input Cards for each Step
      ##----------------------------------------------------

      # use the job_config options to build the Alpgen input card file
      model_input_card_filename,process,njets,inputfilecheck = parse_jobopts(self.job_config,self.evgen_job_opts,job_path)
      
      # read model input card
      model_input_card = AlpgenInputFile()
      model_input_card.read(os.path.join(job_path,model_input_card_filename))

      # filenames for each step
      input_filename_gridgen  = model_input_card.filename_base + '.input.0'
      input_filename_evtgen   = model_input_card.filename_base + '.input.1'
      input_filename_unwght   = model_input_card.filename_base + '.input.2'

      fullproc = process + ' + ' + str(njets) + ' jets'
      gridgen_rate = gridgen_rates[fullproc]['mean']
      if gridgen_rates[fullproc]['n'] != 0:
         gridgen_rate -= gridgen_rates[fullproc]['sigma']
      else:
         gridgen_rate -= gridgen_rate/2.

      # create input for imode 0
      # want about 4 hrs worth of grid gen time so, divide 3/4 iterations and 1/4 full gen
      gridgen_evts = int(4*60*gridgen_rate)
      model_input_card.imode             = 0
      model_input_card.start_with        = 0
      model_input_card.nevt              = int(gridgen_evts*0.75)
      model_input_card.nitr              = 5
      model_input_card.last_nevt         = int(gridgen_evts*0.25/model_input_card.nitr)
      model_input_card.write(os.path.join(job_path,input_filename_gridgen))
      
      # create input for imode 1
      model_input_card.imode             = 1
      model_input_card.start_with        = 2
      model_input_card.nevt              = 0
      model_input_card.nitr              = 0
      model_input_card.last_nevt         = self.nevts
      model_input_card.write(os.path.join(job_path,input_filename_evtgen))
      
      # create input for imode 2
      with open(os.path.join(job_path,input_filename_unwght),'w') as f:
         f.write('2   ! imode\n')
         f.write(model_input_card.filename_base + '   ! base filname\n')
         f.close()


      ##-------------------------------------
      #   Create Grid Generation BalsamJob
      ##-------------------------------------

      #  grid filenames
      grid1 = model_input_card.filename_base + '.grid1'
      grid2 = model_input_card.filename_base + '.grid2'

      gridgen = BalsamJob()
      gridgen.executable         = process + 'gen90_mpi'
      gridgen.executable_args    = input_filename_gridgen
      gridgen.input_files        = [input_filename_gridgen,PDF_FILENAME]
      gridgen.output_files       = [grid1,grid2]
      gridgen.nodes              = 1 # grid generation is serial
      gridgen.processes_per_node = 1 # grid generation is serial
      gridgen.wall_minutes       = 0 # running on a condor cluster with no wall time constraints
      gridgen.username           = os.environ['USER']
      gridgen.target_site        = self.gridgen_site

      ##-------------------------------------
      #   Create Event Generation BalsamJob
      ##-------------------------------------

      preprocess = 'alpgen_presubmit.sh'
      postprocess = 'alpgen_postsubmit.sh'
      output_data_file = model_input_card.filename_base + '.unw.gz'
      output_par_file = model_input_card.filename_base + '_unw.par'

      evgen_rate = evgen_rates[fullproc]
      shower_eff = pythia_shower_eff[njets]

      unw_nevts = self.nevts/shower_eff

      evgen = BalsamJob()
      evgen.executable           = 'alpgenCombo.sh' 
      evgen.executable_args      = process + 'gen90_mpi ' + input_filename_evtgen + ' ' + input_filename_unwght + ' ' + str(16)
      evgen.input_files          = [input_filename_evtgen,
                                    input_filename_unwght,
                                    PDF_FILENAME,
                                    grid1, grid2,
                                   ]
      evgen.output_files         = [
                                    output_data_file,
                                   output_par_file,
                                    'directoryList_before.txt',
                                    'directoryList_after.txt',
                                    postprocess + '.out',
                                    postprocess + '.err',
                                   ]
      evgen.preprocess           = preprocess
      evgen.postprocess          = postprocess
      evgen.postprocess_args     = model_input_card.filename_base
      evgen.nodes                = 1 # grid generation is serial
      evgen.processes_per_node   = 1 # grid generation is serial
      evgen.wall_minutes         = 0 # running on a condor cluster with no wall time constraints
      evgen.username             = os.environ['USER']
      evgen.target_site          = self.evgen_site

      ##-------------------------------------
      #   Create Showering Job: athena_condor_job.py
      ##-------------------------------------

      # ./athena_condor_job.py -g -j 19.2.4.16,AtlasProduction -k /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase -x "--AMITag=e4721 --evgenJobOpts=MC15JobOpts-00-02-22_v0.tar.gz --ecmEnergy=13000 --steering=afterburn --runNumber=361705 --jobConfig=MC15JobOptions/MC15.361705.AlpgenPythiaEvtGen_P2012_ZeeNp5.py --inputGeneratorFile=AlpgenPythia_P2012_ZeeNp5._00001.tar.gz --outputEVNTFile=output_filename.pool.root" -i AlpgenPythia_P2012_ZeeNp5._00001.tar.gz

      ATLASLocalRootBase = '/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase'
      inputGeneratorFile = inputfilecheck + '.tar.gz'
      shower = BalsamJob()
      shower.executable          = 'athena_condor_job.py'
      shower.executable_args     = (' -g -j' + asetup_config + ' -k ' + ATLASLocalRootBase + 
                                    ' -x " --AMITag=' + self.ami_tag + ' --evgenJobOpts=' + self.evgenJobOpts +
                                    ' --ecmEnergy=' + self.ecmEnergy + ' --steering=' + self.steering +
                                    ' --runNumber=' + self.runNumber + ' --jobConfig=' + self.jobConfig +
                                    ' --inputGeneratorFile=' + inputGeneratorFile +
                                    ' -- outputEVNTFile=' + self.output_evnt_file + '"' + 
                                    ' -i ' inputGeneratorFile
                                   )
      shower.preprocess          = 'alpgen_create_athena_input.py'
      shower.preprocess_args     = (' --events=' + output_data_file + 
                                    ' --parameters=' + output_par_file + 
                                    ' --output_base=' + inputfilecheck
                                   )
      shower.nodes               = 1 # shower is serial
      shower.processes_per_node  = 1 # shower is serial
      shower.wall_minutes        = 0 # running on a condor cluster with no wall time constraints
      shower.target_site         = self.shower_site
      shower.input_files         = [
                                    output_data_file,
                                    output_par_file,
                                    
                                    ]




      ##-----------------------
      # create argo job
      ##-----------------------
      argo_job = ArgoJob()
      argo_job.input_url          = 'gsiftp://atlasgridftp02.hep.anl.gov' + job_path
      argo_job.output_url         = 'gsiftp://atlasgridftp02.hep.anl.gov' + job_path
      argo_job.username           = os.environ['USER']
      argo_job.group_identifier   = 'pilot_alpgen_test'
      argo_job.job_status_routing_key = os.environ['HOSTNAME'] + '.' + os.path.basename(job_path)

      argo_job.add_job(gridgen)
      argo_job.add_job(evgen)
      argo_job.add_job(shower)
     


      return argo_job

def parse_jobopts(job_config,evgen_job_opts,working_path):
   # keep track of where we are now
   current_path = os.getcwd()
   # move to working path
   os.chdir(working_path)

   # download the job options tarball
   logger.debug(' downloading ' + EVGEN_JOB_OPTS_URL + evgen_job_opts)
   urllib.urlretrieve(EVGEN_JOB_OPTS_URL + evgen_job_opts,evgen_job_opts)

   # open the tarball
   p = subprocess.Popen(['tar','zxf',evgen_job_opts],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   stdout,stderr = p.communicate()
   if p.returncode != 0:
      logger.error('"tar" returned nonzero value, ' + str(p.returncode) + ', while untarring the evgen_job_opts.')

   # open job config and extract alpgen input card
   card_filename = 'alpgen_tmp_input_card'
   process = ''
   with open(job_config,'rb') as fp:
      job_config_module = imp.load_module('job_config_module',fp,job_config,('.py','rb',imp.PY_SOURCE))
      # write a temporary alpgen input configuration card
      process = job_config_module.alpgen_process
      njets = job_config_module.alpgen_njets
      inputfilecheck = job_config_module.alpgen_inputfilecheck
      with open(card_filename,'w') as outfile:
         outfile.write('''0       ! imode
alpout         ! imode
0              ! start with: 0=new grid, 1=previous warmup grid, 2=previous generation grid
0            0 ! Num events/iteration, Num warm-up iterations
0              ! Num events generated after warm-up
''')
         outfile.write(job_config_module.alpgen_input_card)

   # remove tar ball and what came out of it.
   for thing in os.listdir(os.getcwd()):
      if card_filename in thing:
         pass
      else:
         os.system('rm -rf ' + thing)

   # return to working folder
   os.chdir(current_path)

   return card_filename,process,njets,inputfilecheck





def get_unique_job_path(job_directory = JOB_DIRECTORY):
   uniqueID = str(int(time.time()*1000000))
   while os.path.exists(os.path.join(job_directory,uniqueID)):
      uniqueID = str(int(time.time()*1000000))
   os.mkdir(os.path.join(job_directory,uniqueID))
   return os.path.join(job_directory,uniqueID)

def main():
   
   logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='submit alpgen pilot job to ARGO')
   parser.add_option('-c','--ecmEnergy',dest='ecmEnergy',help='Center of mass energy.',type='int')
   parser.add_option('-r','--runNumber',dest='runNumber',help='Job run number.')
   parser.add_option('-j','--jobConfig',dest='jobConfig',help='A comma-separated list of job configuration script files')
   parser.add_option('-e','--evgenJobOpts',dest='evgenJobOpts',help='Download and install the EvgenJobOpts tarball with the given name.')
   parser.add_option('-o','--outputEVNTFile',dest='outputEVNTFile',help='POOL file into which generated events will be written.')
   parser.add_option('-a','--AMITag',dest='AMITag',help='AMI tag from which this job was defined - this option simply writes the relevant AMI tag value into the output metadata, it does not configure the job.')
   parser.add_option('-s','--steering',dest='steering',help="Steer the transform by manipulating the execution graph before the execution path is calculated. Format is substep:{in,out}{+-}DATA,{in,out}{+-}DATA,... to modify the substep's input/output by adding/removing a data type. e.g. RAWtoESD:in-RDO,in+RDO_TRIG would remove RDO and add RDO_TRIG to the list of valid input datatypes for the RAWtoESD substep..")
   parser.add_option('-g','--gridgen-site',dest='gridgen_site',help='Balsam Site on which to generate the integration grids.')
   parser.add_option('-i','--evgen-site',dest='evgen_site',help='Balsam Site on which to generate the events.')
   parser.add_option('-n','--nevts',dest='nevts',help='Number of events to generate.',type='int')
   parser.add_option('--submit',dest='submit',action='store_true',default=False,help='Must be included to submit the job.')
   options,args = parser.parse_args()


   manditory_args = [
                     'ecmEnergy',
                     'runNumber',
                     'jobConfig',
                     'evgenJobOpts',
                     'outputEVNTFile',
                     'AMITag',
                     'steering',
                     'nevts',
                     'gridgen_site',
                     'evgen_site',
                  ]

   for man in manditory_args:
      if options.__dict__[man] is None:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)

   pilot_alpgen_job = PilotAlpgenArgoJob(options.ecmEnergy,
                            options.runNumber,
                            options.jobConfig,
                            options.evgenJobOpts,
                            options.outputEVNTFile,
                            options.AMITag,
                            options.steering,
                            options.nevts,
                            options.gridgen_site,
                            options.evgen_site,
                           )

   argo_job = pilot_alpgen_job.get_argo_job()
   #logger.info('message = \n' + argo_job.serialize())


   if options.submit:
      # load message interface, but choose the development version if we are developing:
      sys.path.append('/users/hpcusers/argobalsam/production/argobalsam/common_core')
      from MessageInterface import MessageInterface
      
      mi = MessageInterface()
      mi.host = 'atlasgridftp02.hep.anl.gov'
      mi.port = 5671
      mi.ssl_cert = os.environ['X509_USER_CERT'] #'/users/hpcusers/balsam/gridsecurity/jchilders/xrootdsrv-cert.pem'
      mi.ssl_key  = os.environ['X509_USER_KEY'] #'/users/hpcusers/balsam/gridsecurity/jchilders/xrootdsrv-key.pem'
      mi.ssl_ca_certs = os.environ['X509_CA_CERTS'] #'/users/hpcusers/balsam/gridsecurity/jchilders/cacerts.pem'
      mi.exchange_name = 'argo_users'
      
      logger.info('opening connection')
      mi.open_blocking_connection()
      routing_key = 'argo_job'
      logger.info(' sending msg ')
      mi.send_msg(argo_job.serialize(),routing_key)
      logger.info(' done sending')
      mi.close()
      logger.info(' closing connection')
   else:
      logger.info(' not submitting job ')

   '''
   if options.enable_status_queue:
      argo_job.job_status_routing_key = 'test_job_status' #'status_' + jobID
      logger.info('setting up job status queue: ' + argo_job.job_status_routing_key)
      mi = MessageInterface()
      mi.host = 'atlasgridftp02.hep.anl.gov'
      mi.port = 5671
      mi.ssl_cert = os.environ['X509_USER_CERT'] #'/users/hpcusers/balsam/gridsecurity/jchilders/xrootdsrv-cert.pem'
      mi.ssl_key  = os.environ['X509_USER_KEY'] #'/users/hpcusers/balsam/gridsecurity/jchilders/xrootdsrv-key.pem'
      mi.ssl_ca_certs = os.environ['X509_CA_CERTS'] #'/users/hpcusers/balsam/gridsecurity/jchilders/cacerts.pem'
      mi.exchange_name = 'argo_users'
      if options.dev:
         mi.exchange_name = 'argo_users_dev'
      mi.open_blocking_connection()
      mi.create_queue(argo_job.job_status_routing_key,argo_job.job_status_routing_key)
      mi.close()
   '''






if __name__ == '__main__':
   main()


