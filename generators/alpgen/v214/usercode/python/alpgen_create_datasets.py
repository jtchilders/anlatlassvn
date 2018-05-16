#!/usr/bin/env python
import os,sys,optparse,logging

import AlpgenUnwFile_v2
import AlpgenUnwParFile

# ---- PYTHIA Default settings
DEFAULT_PYTHIA_BINARY = '/users/hpcusers/svn/generators/pythia/v8180/usercode/alpToHepmc/trunk/runPythiaOnAlpgen'
DEFAULT_PYTHIA_XMLDOC = '/users/hpcusers/svn/generators/pythia/v8180/unedited/xmldoc'
DEFAULT_PYTHIA_JET_MATCHING_INCLUSIVE = False
DEFAULT_PYTHIA_NEVTS = 5000
DEFAULT_PYTHIA_OUTPUT_FILENAME = 'pythia_output.hepmc.tmp'

DEFAULT_UNW_EVTS_PER_FILE = -1
DEFAULT_SHWR_EVTS_PER_FILE = 5700
DEFAULT_FILES_PER_SUBSET = 1000

DEFAULT_READ_WRITE_SIZE = 1024*1024

DEFAULT_NFILES       = -1
DEFAULT_FIRST_EVENT  = 0
DEFAULT_FIRST_FILE   = 0

def main():
   logging.basicConfig(level=logging.INF,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='Use to reorganize Alpgen Unweighted events and parameter files. Clearly the inputs should all be of the same type. This script does not check for mixing samples.')

   # input/output
   parser.add_option('-i','--input',dest='input',help='This provides the input Alpgen Unweighted basenames. For example, providing "/some/path/to/alpout" assumes that the unweighted events are in "/some/path/to/alpout.unw" and the associated parameters file is in "/some/path/to/alpout_unw.par". This can be a string of comma seperated basenames.')
   parser.add_option('-o','--output-base',dest='output_basename',help='The basename to use for the output files. For example, if "/some/path/to/output_" is given the output files will be "/some/path/to/output_<file-number>.unw" and "/some/path/to/output_<file-number>_unw.par".')

   # numbers per ...
   parser.add_option('-a','--unw-evts-per-file',dest='unw_evts_per_file',help='The number of unweighted events to include in each file. If not specified, this number will be calculated from running Pythia and extracting the showering effeciency. Then calculating the number of events need such that Pythia outputs <shr-evts-per-file>. [default='+str(DEFAULT_UNW_EVTS_PER_FILE)+']',type='int',default=DEFAULT_UNW_EVTS_PER_FILE)
   parser.add_option('-b','--shwr-evts-per-file',dest='shwr_evts_per_file',help='The target number of showered events per dataset file. [default='+str(DEFAULT_SHWR_EVTS_PER_FILE) + ']',type='int',default=DEFAULT_SHWR_EVTS_PER_FILE)
   parser.add_option('-c','--files-per-subset',dest='files_per_subset',help='The number of files to put in each sub dataset. [default='+str(DEFAULT_FILES_PER_SUBSET)+']',type='int',default=DEFAULT_FILES_PER_SUBSET)
   parser.add_option('-n','--nevts-shwr-final',dest='nevts_shwr_final',help='The total number of showered events to when completed.')

   # dataset name
   parser.add_option('-d','--dataset-name',dest='dataset_name',help='This is the base name for the dataset which will be used for the file names inside the datasets and the container name. The sub-dataset names will be derived from this name using "<dataset-name>_iXY" where X is the upload attempt and Y is the sub-datset number.')


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
   

   options,args = parser.parse_args()

   
   manditory_args = [
                     'input',
                     'output_basename',
                     'unw_evts_per_file',
                     'shwr_evts_per_file',
                     'files_per_subset',
                     'nevts_shwr_final',
                     'dataset_name',
                     'pythia_binary',
                     'pythia_inclusive',
                     'pythia_nevts',
                     'pythia_xmldoc_path',
                    ]

   for man in manditory_args:
      if not options.__dict__[man]:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         return -100

   input_basenames = input.split(',')

   alpgen_create_datasets(
                           input_basenames      = input_basenames,
                           output_basename      = options.output_basename,

                           unw_evts_per_file    = options.unw_evts_per_file,
                           shwr_evts_per_file   = options.shwr_evts_per_file,
                           files_per_subset     = options.files_per_subset,
                           nevts_shwr_final     = options.nevts_shwr_final,

                           dataset_name         = options.dataset_name,

                           pythia_binary        = options.pythia_binary,
                           pythia_inclusive     = options.pythia_inclusive,
                           pythia_nevts         = options.pythia_nevts,
                           pythia_xmldoc_path   = options.pythia_xmldoc_path,
                        )

def alpgen_create_datasets(
                           input_basenames,
                           output_basename,

                           dataset_name,
                           nevts_shwr_final,

                           unw_evts_per_file    = DEFAULT_UNW_EVTS_PER_FILE,
                           shwr_evts_per_file   = DEFAULT_SHWR_EVTS_PER_FILE,
                           files_per_subset     = DEFAULT_FILES_PER_SUBSET,

                           pythia_binary        = DEFAULT_PYTHIA_BINARY,
                           pythia_inclusive     = DEFAULT_PYTHIA_JET_MATCHING_INCLUSIVE,
                           pythia_nevts         = DEFAULT_PYTHIA_NEVTS,
                           pythia_xmldoc_path   = DEFAULT_PYTHIA_XMLDOC,
                           ):
   
   if len(input_basenames) <= 0:
      logger.error('No inputs provided')
      raise Exception('No Input Files provided.')

   # record current working directory
   CWD = os.cwd()
   # create working directory
   os.mkdir(dataset_name)
   # move to new directory
   os.chdir(dataset_name)

   # get the number of events to include per file
   if unw_evts_per_file == DEFAULT_UNW_EVTS_PER_FILE:
      input_filename = input_basenames[0] + '.unw'
      unw_evts_per_file = get_unw_evts_per_file(input_filename,
                              pythia_binary,pythia_inclusive,
                              pythia_nevts,pythia_xmldoc_path)


   # reorganize the input files into sub-datasets
   dataset_counter      = 0
   total_file_counter   = 0


   output_file_counter  = 0
   output_unw_file      = None
   output_unw_filename  = ''
   output_par_file      = None
   output_par_filename  = ''

   for basename in input_basenames:
      input_unw_filename = basename + '.unw'
      input_par_filename = basename + '_unw.par'

      if not os.path.exists(par_filename):
         raise Exception('Unweighted par file does not exists: ' + par_filename)
      input_par_file = AlpgenUnwParFile.AlpgenUnwParFile.read_file(par_filename)

      if not os.path.exists(unw_filename):
         raise Exception('Unweighted events file does not exists: ' + unw_filename)
      input_unw_file = AlpgenUnwFile_v2.AlpgenUnwFile(unw_filename,AlpgenParameters.event_sizes[])

      # need to create new output files if they aren't set
      if ouptut_unw_file is None and output_par_file is None: 
         output_unw_filename = get_output_filename(output_basename,output_file_counter,'.unw')
         output_par_filename = get_output_filename(output_basename,output_file_counter,'_unw.par')

         bytes_per_event = input_par_file.get_bytes_per_event()
         output_unw_file = AlpgenUnwFile_v2.AlpgenUnwFile(output_unw_filename,bytes_per_event)

      # copy input unw file to output unw file in chunks



def get_output_filename(basename,counter,extension):
   return '%s%08i%s' %(basename,counter,extension)


def get_unw_evts_per_file(
                           input_filename,
                           binary = DEFAULT_PYTHIA_BINARY,
                           inclusive = DEFAULT_PYTHIA_JET_MATCHING_INCLUSIVE,
                           shower_nevts = DEFAULT_PYTHIA_NEVTS,
                           xmldoc_path = DEFAULT_PYTHIA_XMLDOC,
                           output_filename = DEFAULT_PYTHIA_OUTPUT_FILENAME,
                         ):
   cmd = pythia_binary
   if pythia_inclusive:
      cmd += ' -i' # turn on inclusive jets
   cmd += ' -n ' + str(shower_nevts) # set number of events to process
   cmd += ' -x ' + xmldoc_path # set xmldoc path
   cmd += ' -o ' + output_filename # set output filename (to be deleted later)
   cmd += ' -a ' + input_filename # add input file basename
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


   p.wait()
   if p.returncode != 0:
      logger.error('pythia exited with non-zero return code: ' + str(p.returncode) + '. STDERR = ')
      for line in p.stderr:
         logger.error(line[0:-1])
      if os.path.exists(output_filename):
         os.remove(output_filename)
      raise Exception('Pythia Exited with none-zero return code')

   # remove pythia file
   os.remove(output_filename)

   # calculate the number of unw events per file rounded to the nearest 1000
   unw_evt_per_file = ( int( num_events / pythia_efficiency / 1000. ) + 1) * 1000
   logger.info(' unweighted events per file:  %i' % unw_evt_per_file)
   logger.info(' pythia efficiency: %10.8e' % pythia_efficiency)

   return unw_evt_per_file


if __name__ == "__main__":
   sys.exit(main())