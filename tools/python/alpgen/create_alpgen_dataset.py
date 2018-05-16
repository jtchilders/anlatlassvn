#!/usr/bin/env python
import os,sys,optparse,logging,subprocess
logger = logging.getLogger(__name__)
logger.name = os.path.basename(sys.argv[0])

sys.path.append('/users/hpcusers/svn/generators/alpgen/v214/usercode/python')
import split_alpgen_unw
from AlpgenUnwParFile import AlpgenUnwParFile

TMP_WORKING_FOLDER = 'tmp_create_alpgen_dataset'

DEFAULT_ALPGEN_BASE = os.path.join(os.environ['PWD'],'alpout')
DEFAULT_NUM_EVENTS = 6000
DEFAULT_PRESHOWER_EVENTS = -1
DEFAULT_PYTHIA_BINARY = '/users/hpcusers/svn/generators/pythia/v8180/usercode/alpToHepmc/trunk/runPythiaOnAlpgen'
DEFAULT_PYTHIA_XMLDOC = '/users/hpcusers/svn/generators/pythia/v8180/unedited/xmldoc'
DEFAULT_PYTHIA_JET_MATCHING_INCLUSIVE = False
DEFAULT_PYTHIA_NEVTS = 5000
DEFAULT_PYTHIA_OUTPUT_FILENAME = 'pythia_output.hepmc.tmp'


DEFAULT_SPLIT_MAX_FILES = split_alpgen_unw.DEFAULT_MAXFILES
DEFAULT_SPLIT_OUTPUT_OFFSET = split_alpgen_unw.DEFAULT_OUTPUT_OFFSET
DEFAULT_SPLIT_MAKE_TARBALLS = True
DEFAULT_SPLIT_SKIP_LAST = True
DEFAULT_SPLIT_TARBALL_BASE = split_alpgen_unw.ALPGEN_ATHENA_FILEBASE
DEFAULT_SPLIT_FIRST_EVENT = split_alpgen_unw.DEFAULT_FIRST_EVENT
DEFAULT_SPLIT_FIRST_FILE = split_alpgen_unw.DEFAULT_FIRST_FILE

ASSUMED_SHOWERED_EVENTS_PER_POOL_FILE = 5000

PYTHIA_SEARCH_STRING = '| sum                                                |'
RUCIO_SCRIPT_FILENAME   = 'rucio_upload.sh'
RUCIO_UPLOAD_SITE     = 'ANLASC_SCRATCHDISK'



def main():
   logging.basicConfig(level=logging.INFO)

   parser = optparse.OptionParser(description=' Create and Upload a datset of alpgen unweighted events (and par) files for ATLAS consumption. The input files should all be of the same type, do not mix processes. Pythia is used to characterize the jet matching efficiency to guide the number of events to include in each output file.')
   parser.add_option('-i','--input',dest='input',
      help='This should a be a comma separated list of alpgen base file names, [default = "' + DEFAULT_ALPGEN_BASE + '"]',
      default=DEFAULT_ALPGEN_BASE)
   parser.add_option('-n','--num-events',dest='num_events',
      help='The number of events per file after showering. This is used to calculate how many unweighted events should be included in each file that is put into the dataset. [default = ' + str(DEFAULT_NUM_EVENTS) + ']',
      default=DEFAULT_NUM_EVENTS,type='int')
   parser.add_option('-e','--events',dest='events',help='Number of events needed in the dataset after showering.',type='int')
   parser.add_option('--num-preshower-evts',dest='num_preshower_evts',help='If set, overrides running Pythia to find the number of unweighted events needed per file and uses the passed value',type='int',default=DEFAULT_PRESHOWER_EVENTS)
   parser.add_option('--retry-upload',dest='retry_upload',help='Only retry uploading to grid',action='store_true',default=False)
   parser.add_option('--production-role',dest='production_role',help='Use RUCIO_ACCOUNT=phys-gener instead of "childers"',action='store_true',default=False)
   
   # ------------------------------------
   # pythia related options
   parser.add_option('--pythia-binary',dest='pythia_binary',
      help='This script runs pythia over 5000 events to shower them and get an efficiency of parton-jet matching for the output file. That way this efficiency is used to calculate the number of unweighted events to put in the pre-showered event files. [default = "' + DEFAULT_PYTHIA_BINARY + '"]',
      default=DEFAULT_PYTHIA_BINARY)
   parser.add_option('--pythia-inclusive',dest='pythia_inclusive',
      help='Pythia option to turn on inclusive jet matching.[default='+str(DEFAULT_PYTHIA_JET_MATCHING_INCLUSIVE)+']',
      default=DEFAULT_PYTHIA_JET_MATCHING_INCLUSIVE,action='store_true')
   parser.add_option('--pythia-nevts',dest='pythia_nevts',type='int',
      help='Pythia option: how many events to shower. [default=' + str(DEFAULT_PYTHIA_NEVTS) + ']',
      default = DEFAULT_PYTHIA_NEVTS)
   parser.add_option('--pythia-xmldoc',dest='pythia_xmldoc_path',
      help='Path to pythia "xmldoc" area. [default="' + DEFAULT_PYTHIA_XMLDOC + '"]',
      default = DEFAULT_PYTHIA_XMLDOC)
   
   # ------------------------------------
   # alpgen split unw file options
   parser.add_option('-d','--dataset-name',dest='dataset_name',
      help='Name of the dataset to be uploaded. The unweighted event files will also contain this string to ensure uniqueness.')

   parser.add_option('--split-max-files',dest='split_max_files',
      help='The maximum number of files to produce when splitting the alpgen unweighted files.[default = ' + str(DEFAULT_SPLIT_MAX_FILES) + ']',
      default = DEFAULT_SPLIT_MAX_FILES)
   parser.add_option('--split-output-offset',dest='split_output_offset',
      help='Splitting Alpgen files option: The offset to start with when creating numbered outputfiles. [DEFAULT=' + str(DEFAULT_SPLIT_OUTPUT_OFFSET)+']',
      type='int',default=DEFAULT_SPLIT_OUTPUT_OFFSET)
   parser.add_option('--split-make-tarballs',dest='split_make_tarballs',help='Splitting Alpgen option: make the tarballs, and for the moment this only makes the tarballs. [default=' + str(DEFAULT_SPLIT_MAKE_TARBALLS) + ']',
      action='store_true',default=DEFAULT_SPLIT_MAKE_TARBALLS)
   parser.add_option('--split-skip-last',dest='split_skip_last',
      help='Splitting Alpgen option: skip the last file when splitting as it will have fewer events than wanted. [default = ' + str(DEFAULT_SPLIT_SKIP_LAST) + ']',
      action='store_false',default=DEFAULT_SPLIT_SKIP_LAST)
   parser.add_option('--split-tarball-base',dest='split_tarball_base',
      help='Splitting Alpgen option: tarball name base. Defaults to dataset name.')
   parser.add_option('--split-first-event',dest='split_first_event',
      help='Splitting Alpgen option: first event to start with. [default=' + str(DEFAULT_SPLIT_FIRST_EVENT) + ']',
      type='int',default=DEFAULT_SPLIT_FIRST_EVENT)
   parser.add_option('--split-first-file',dest='split_first_file',
      help='Splitting Alpgen option: first file to start with. [default=' + str(DEFAULT_SPLIT_FIRST_FILE) + ']',
      type='int',default=DEFAULT_SPLIT_FIRST_FILE)

   parser.add_option('--no-upload',dest='dataset_upload',help='For testing, skip dataset uploading.',action='store_false',default=True)

   options,args = parser.parse_args()

   if options.input is None and not options.retry_upload:
      parser.error('Must specify -i')
   elif options.events is None and not options.retry_upload:
      parser.error('Must specify -e')
   elif options.dataset_name is None:
      parser.error('Must specify -d')

   if options.split_tarball_base is None:
      options.split_tarball_base = options.dataset_name

   return create_alpgen_dataset(
                                input                 = options.input,
                                events                = options.events,
                                dataset_name          = options.dataset_name,
                                num_events            = options.num_events,
                                num_preshower_evts    = options.num_preshower_evts,
                                pythia_binary         = options.pythia_binary,
                                pythia_inclusive      = options.pythia_inclusive,
                                pythia_nevts          = options.pythia_nevts,
                                pythia_xmldoc_path    = options.pythia_xmldoc_path,
                                split_max_files       = options.split_max_files,
                                split_output_offset   = options.split_output_offset,
                                split_make_tarballs   = options.split_make_tarballs,
                                split_skip_last       = options.split_skip_last,
                                split_tarball_base    = options.split_tarball_base,
                                split_first_event     = options.split_first_event,
                                split_first_file      = options.split_first_file,
                                dataset_upload        = options.dataset_upload,
                                retry_upload          = options.retry_upload,
                                production_role       = options.production_role,
                               )

def create_alpgen_dataset(
                          input,
                          events,
                          dataset_name,
                          num_events = DEFAULT_NUM_EVENTS,
                          num_preshower_evts = DEFAULT_PRESHOWER_EVENTS,
                          pythia_binary = DEFAULT_PYTHIA_BINARY,
                          pythia_inclusive = DEFAULT_PYTHIA_JET_MATCHING_INCLUSIVE,
                          pythia_nevts = DEFAULT_PYTHIA_NEVTS,
                          pythia_xmldoc_path = DEFAULT_PYTHIA_XMLDOC,
                          split_max_files = DEFAULT_SPLIT_MAX_FILES,
                          split_output_offset = DEFAULT_SPLIT_OUTPUT_OFFSET,
                          split_make_tarballs = DEFAULT_SPLIT_MAKE_TARBALLS,
                          split_skip_last = DEFAULT_SPLIT_SKIP_LAST,
                          split_tarball_base = DEFAULT_SPLIT_TARBALL_BASE,
                          split_first_event = DEFAULT_SPLIT_FIRST_EVENT,
                          split_first_file = DEFAULT_SPLIT_FIRST_FILE,
                          dataset_upload = True,
                          retry_upload = False,
                          production_role = False,
                         ):
   logger.info(' input basenames: ' + input)
   logger.info(' target number of events per file after showering: ' + str(num_events))
   logger.info(' number of final showered events in the output dataset: ' + str(events))
   logger.info(' number of preshowered events requested: ' + str(num_preshower_evts))
   logger.info(' dataset name: ' + dataset_name)
   logger.info(' pythia binary: ' + pythia_binary)
   jet_matching = None
   if pythia_inclusive:
      logger.info('     pythia jet-matching is inclusive')
      jet_matching = AlpgenUnwParFile.INCLUSIVE_JET_MATCHING
   else:
      logger.info('     pythia jet-matching is exclusive')
      jet_matching = AlpgenUnwParFile.EXCLUSIVE_JET_MATCHING
   logger.info('     pythia nevts: ' + str(pythia_nevts))

   # move to working directory
   initial_working_directory = os.environ['PWD']
   if not os.path.exists(TMP_WORKING_FOLDER):
      os.mkdir(TMP_WORKING_FOLDER)
   os.chdir(TMP_WORKING_FOLDER)

   # create a list of input alpgen basenames to loop over
   input_alpgen_basenames = input.replace(' ','').split(',')
   if len(input_alpgen_basenames) == 0:
      parser.error('Must specify input alpgen files')

   # ---------------------------------------------------
   #  RUN PYTHIA to get the showering efficiency
   #   only need to run this on the first file
   # ---------------------------------------------------

   unw_evt_per_file = num_preshower_evts
   if unw_evt_per_file == DEFAULT_PRESHOWER_EVENTS and not retry_upload:

      cmd = pythia_binary
      if pythia_inclusive:
         cmd += ' -i' # turn on inclusive jets
      cmd += ' -n ' + str(pythia_nevts) # set number of events to process
      cmd += ' -x ' + pythia_xmldoc_path # set xmldoc path
      cmd += ' -o ' + DEFAULT_PYTHIA_OUTPUT_FILENAME # set output filename (to be deleted later)
      cmd += ' -a ' + input_alpgen_basenames[0] # add input file basename
      p = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      pythia_tried_events = 0
      pythia_selected_events = 0
      pythia_accepted_events = 0
      pythia_efficiency = 0.
      for line in p.stdout:
         logger.info(line[0:-1])
         if PYTHIA_SEARCH_STRING in line:
            split = line.split()
            #logger.info(str(split))
            pythia_tried_events = int(split[3])
            pythia_selected_events = int(split[4])
            pythia_accepted_events = int(split[5])
            if pythia_accepted_events != pythia_nevts:
               logger.warning(' Number of pythia events accepted = ' + str(pythia_accepted_events) + ' when ' + str(pythia_nevts) + ' were requested. Pythia tried ' + str(pythia_tried_events) + '.')
            pythia_efficiency = 1.*pythia_accepted_events/pythia_tried_events
            logger.info(' Pythia Showering Efficiency = ' + str(pythia_efficiency))


      p.wait()
      if p.returncode != 0:
         logger.error('pythia exited with non-zero return code: ' + str(p.returncode) + '. STDERR = ')
         for line in p.stderr:
            logger.error(line[0:-1])
         return -1

      # remove pythia file
      os.remove(DEFAULT_PYTHIA_OUTPUT_FILENAME)

      # calculate the number of unw events per file rounded to the nearest 1000
      unw_evt_per_file = ( int( num_events / pythia_efficiency / 1000. ) + 1) * 1000
   logger.info(' spliting unweighted events into ' + str(unw_evt_per_file) + ' events per file.')

   # ---------------------------------------------------
   # SPLIT UNW FILES
   # ---------------------------------------------------

   
   if not retry_upload:
      # calculate the number of files needed
      updated_split_max_files = int( 1. * events / ASSUMED_SHOWERED_EVENTS_PER_POOL_FILE ) + 1
      # if the max files command line options has been set away from its default, then override
      # this calculation
      if split_max_files is not DEFAULT_SPLIT_MAX_FILES:
         updated_split_max_files = split_max_files

      # split files
      split_alpgen_unw.TMP_FOLDER = dataset_name
      ret = split_alpgen_unw.split_alpgen_unw(
                                          input_alpgen_basenames,
                                          unw_evt_per_file,
                                          dataset_name + '.',
                                          updated_split_max_files,
                                          split_output_offset,
                                          split_make_tarballs,
                                          split_skip_last,
                                          split_tarball_base,
                                          split_first_event,
                                          split_first_file,
                                          dataset_name,
                                          tarball_only = split_make_tarballs,
                                          update_par = True,
                                          jet_matching = jet_matching,
                                        )

      if ret < 0:
         logger.error(' Error in split_alpgen_unw. Returned = ' + str(ret))
         return ret


      # if the number of output files was not enough, exit with an error
      if ret < updated_split_max_files:
         logger.error('Not enough events to produce the needed dataset. Needed number of files: ' + str(updated_split_max_files) + ', produced: ' + str(ret) )
         return -1

   # ---------------------------------------------------
   # UPLOAD DATASET
   # ---------------------------------------------------

   # create rucio script
   write_rucio_script(production_role)

   # create the script input
   dot_split = dataset_name.split('.')
   rucio_dataset_scope = dot_split[0] + '.' + dot_split[1]
   dataset_with_scope = rucio_dataset_scope + ':' + dataset_name

   # run script
   cmd = 'sh ' + RUCIO_SCRIPT_FILENAME + ' ' + RUCIO_UPLOAD_SITE + ' ' + rucio_dataset_scope + ' ' + dataset_with_scope + ' ' + dataset_name + '/'
   logger.info('running command: ' + cmd)
   
   if dataset_upload:
      for i in range(5): # retry up to five times
         p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
         logger.info('STDOUT from subprocess:')
         for line in p.stdout:
            logger.info(line[0:-1])
         p.wait()
         if p.returncode != 0:
            logger.error(' script returned a non-zero value, ' + str(p.returncode) + '; STDERR = ')
            for line in p.stderr:
               logger.error(line[0:-1])
            if p.returncode == 1: # seen with server errors
               continue # try again
            elif p.returncode == 2: # seen when command line syntax is wrong
               break # need to fix syntax
            elif p.returncode < 0: # -15 seen when file size/checksum are off and when upload just hangs. -9 for 'bad gateway'
               continue # try again
            return -1
         else:
            logger.info('STDERR from subprocess:')
            for line in p.stderr:
               logger.info(line[0:-1])
         break # don't try again if all goes well
   
   logger.info(' uploaded dataset: ' + dataset_name)
   if not retry_upload:
      logger.info(' events per file: ' + str(unw_evt_per_file))
      logger.info(' number of files: ' + str(updated_split_max_files))

   # return to original folder
   os.chdir(initial_working_directory)
   return 0
   

def write_rucio_script(
                        production_role         = False,
                        atlas_local_root_base   = '/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase',
                        rucio_script_filename   = RUCIO_SCRIPT_FILENAME,
                        rucio_account_user      = 'childers',
                        rucio_account_prod      = 'phys-gener',
                      ):
   
   rucio_accout = rucio_account_user
   if production_role:
      rucio_accout = rucio_account_prod

   rucio_script_content = """#!/usr/bin/env bash
RESOURCE=$1
SCOPE=$2
DATASET_NAME=$3
SOURCE_DIRECTORY=$4
# setup ASetup environment
export ATLAS_LOCAL_ROOT_BASE=""" + atlas_local_root_base + """
# setup ASetup environment
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
echo ' setting up dq2 '
export RUCIO_ACCOUNT=""" + rucio_accout + """
localSetupRucioClients
echo ' voms proxy status: '
voms-proxy-info -all
echo ' running rucio upload with: '
echo '   RESOURCE          = ' $RESOURCE
echo '   SCOPE             = ' $SCOPE
echo '   DATASET_NAME      = ' $DATASET_NAME
echo '   SOURCE_DIRECTORY  = ' $SOURCE_DIRECTORY
echo '   RUCIO_ACCOUNT     = ' $RUCIO_ACCOUNT
rucio -v upload --rse $RESOURCE --scope $SCOPE $DATASET_NAME $SOURCE_DIRECTORY 2>&1
echo ' rucio exit code:' $?
"""
   rucio_script = open(rucio_script_filename,'w')
   rucio_script.write(rucio_script_content)
   rucio_script.close()

if __name__ == "__main__":
   sys.exit(main())
