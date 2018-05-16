#!/usr/bin/env python
import os,sys,optparse,logging,glob,copy,subprocess,shutil
logger = logging.getLogger(__name__)
import AlpgenUnwFile_v2
import AlpgenUnwParFile
import AlpgenParameters
sys.path.append('/users/hpcusers/svn/tools/python')
import rucio_bulk_upload

MAX_SUBSETS = 9 # this comes from the 'mc15_v1_iXY' notation. Y = the number of subsets
MAX_FILES_PER_SUBSET = 99999 # this comes from the 'mc15_v1._00000' notation. only 5 digits in the file counter

DEFAULT_NFILES             = -1
DEFAULT_FIRST_EVENT        = 0
DEFAULT_FIRST_FILE         = 0
DEFAULT_SUBSETS            = 1
DEFAULT_COLLECTION_NAME    = 'alpgen_data_collection'

DEFAULT_UNW_EVTS_PER_FILE = -1
DEFAULT_SHWR_EVTS_PER_FILE = 5700
DEFAULT_FILES_PER_SUBSET = MAX_FILES_PER_SUBSET

DEFAULT_SUBSET_INDEX       = -9

# this sets the scale of the read/write block size
# it is rounded up based on the event size, that way
# an integer number of events are read/written at one time
DEFAULT_READ_WRITE_BLOCK_SIZE = 1024*1024 

DEFAULT_OUTPUT_FILE_INDEX_OFFSET = 0



# ---- PYTHIA Default settings
USE_PYTHIA_6 = True
DEFAULT_PYTHIA_BINARY_6 = '/users/hpcusers/svn/generators/pythia/v6428/usercode/alpToHepmc/alpToHepmc'
DEFAULT_PYTHIA_BINARY_8 = '/users/hpcusers/svn/generators/pythia/v8180/usercode/alpToHepmc/trunk/runPythiaOnAlpgen'
if USE_PYTHIA_6:
   DEFAULT_PYTHIA_BINARY = DEFAULT_PYTHIA_BINARY_6
else:
   DEFAULT_PYTHIA_BINARY = DEFAULT_PYTHIA_BINARY_8
DEFAULT_PYTHIA_XMLDOC = '/users/hpcusers/svn/generators/pythia/v8180/unedited/xmldoc'
DEFAULT_PYTHIA_JET_MATCHING_INCLUSIVE = False
DEFAULT_PYTHIA_NEVTS = 5000
DEFAULT_PYTHIA_OUTPUT_FILENAME = '/tmp/pythia_output.hepmc.tmp.pid' + str(os.getpid())

def main():
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='Use to reorganize Alpgen Unweighted events and parameter files. Clearly the inputs should all be of the same type. This script does not check for mixing samples.')

   parser.add_option('-i','--input',dest='input',help='This provides the input Alpgen Unweighted basenames. For example, providing "/some/path/to/alpout" assumes that the unweighted events are in "/some/path/to/alpout.unw" and the associated parameters file is in "/some/path/to/alpout_unw.par". This can be a string of comma seperated basenames. This will also look for "/some/path/to/alpout.unw.gz" if the ".unw" does not exist.')
   #parser.add_option('-o','--output-path',dest='output_path',help='The path to use for the output files.')

   parser.add_option('-c','--collection-name',dest='collection_name',help='The folder name where the subsets will be created. [default='+DEFAULT_COLLECTION_NAME+']',default=DEFAULT_COLLECTION_NAME)
   #parser.add_option('-d','--dataset-name',dest='dataset_name',help='Subsets will be named <dataset-name>. If not speficied the <collection-name>_iXY will be used.')

   parser.add_option('--subset-index',dest='subset_index',help='If changed from the default of '+str(DEFAULT_SUBSET_INDEX)+' it will override the automatic detection of previous subsets.',default=DEFAULT_SUBSET_INDEX,type='int')

   parser.add_option('-a','--first-event',dest='first_event',help='The first event to start with. This defaults to the first event of the first file. [default='+str(DEFAULT_FIRST_EVENT)+']',type='int',default=DEFAULT_FIRST_EVENT)
   parser.add_option('-b','--first-file',dest='first_file',help='The first file to start with. This defaults to the first  file. [default='+str(DEFAULT_FIRST_FILE)+']',type='int',default=DEFAULT_FIRST_FILE)

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
   
   # numbers per ...
   parser.add_option('-d','--unw-evts-per-file',dest='unw_evts_per_file',help='The number of unweighted events to include in each file. If not specified, this number will be calculated from running Pythia and extracting the showering effeciency. Then calculating the number of events need such that Pythia outputs <shr-evts-per-file>. [default='+str(DEFAULT_UNW_EVTS_PER_FILE)+']',type='int',default=DEFAULT_UNW_EVTS_PER_FILE)
   parser.add_option('-e','--shwr-evts-per-file',dest='shwr_evts_per_file',help='The target number of showered events per dataset file. [default='+str(DEFAULT_SHWR_EVTS_PER_FILE) + ']',type='int',default=DEFAULT_SHWR_EVTS_PER_FILE)
   parser.add_option('-g','--files-per-subset',dest='files_per_subset',help='The number of files to put in each sub dataset. [default='+str(DEFAULT_FILES_PER_SUBSET)+']',type='int',default=DEFAULT_FILES_PER_SUBSET)
   parser.add_option('-n','--nevts-shwr-final',dest='nevts_shwr_final',help='The total number of showered events to when completed.',type='int')

   # additional controls
   parser.add_option('','--output-file-offset',dest='output_file_index_offset',help='Output files are numbered starting from 0 by default, but this allows the user to change the starting file number.[default=' + str(DEFAULT_OUTPUT_FILE_INDEX_OFFSET) +']',type='int',default=DEFAULT_OUTPUT_FILE_INDEX_OFFSET)

   # for uploading
   parser.add_option('','--upload2grid',dest='upload',action='store_true',default=False,help='This will upload the dataset to the grid, but the proxy must exist and be valid')
   parser.add_option('','--proxy',dest='proxy',help='the location of the proxy certificate')
   parser.add_option('','--db-entries',dest='db_entries',default='',help='a comma separated list of db entry ids to update when the uploading begins and when it ends so as to keep track of entries that were uploaded.')

   options,args = parser.parse_args()

   manditory_args = [
                     'input',
                     #'output_basename',

                     'collection_name',
                     'subset_index',

                     'pythia_binary',
                     'pythia_inclusive',
                     'pythia_nevts',
                     'pythia_xmldoc_path',

                     'unw_evts_per_file',
                     'shwr_evts_per_file',
                     'files_per_subset',
                     'nevts_shwr_final',

                     'output_file_index_offset',
                    ]
   
   for man in manditory_args:
      if options.__dict__[man] is None:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         return -100


   if options.first_event != DEFAULT_FIRST_EVENT and options.first_file != DEFAULT_FIRST_FILE:
      logger.error('Cannot specify both the first event and the first file. Not yet supported.')
      return -100

   input_basenames = options.input.split(',')

   alpgen_reorganize_unw(
                           input_basenames      = input_basenames,
                           collection_name      = options.collection_name,
                           nevts_shwr_final     = options.nevts_shwr_final,

                           first_event          = options.first_event,
                           first_file           = options.first_file,

                           pythia_binary        = options.pythia_binary,
                           pythia_inclusive     = options.pythia_inclusive,
                           pythia_nevts         = options.pythia_nevts,
                           pythia_xmldoc_path   = options.pythia_xmldoc_path,

                           unw_evts_per_file    = options.unw_evts_per_file,
                           shwr_evts_per_file   = options.shwr_evts_per_file,
                           files_per_subset     = options.files_per_subset,

                           output_file_index_offset = options.output_file_index_offset,

                           subset_index         = options.subset_index,

                           upload               = options.upload,
                           proxy                = options.proxy,
                           db_entries           = options.db_entries
                        )

class MaxNumSubsetExceeded(Exception): pass
class MaxFilesPerSubsetExceeded(Exception): pass
def alpgen_reorganize_unw(
                           input_basenames,
                           collection_name,
                           nevts_shwr_final,
                           
                           first_event          = DEFAULT_FIRST_EVENT,
                           first_file           = DEFAULT_FIRST_FILE,

                           pythia_binary        = DEFAULT_PYTHIA_BINARY,
                           pythia_inclusive     = DEFAULT_PYTHIA_JET_MATCHING_INCLUSIVE,
                           pythia_nevts         = DEFAULT_PYTHIA_NEVTS,
                           pythia_xmldoc_path   = DEFAULT_PYTHIA_XMLDOC,

                           unw_evts_per_file    = DEFAULT_UNW_EVTS_PER_FILE,
                           shwr_evts_per_file   = DEFAULT_SHWR_EVTS_PER_FILE,
                           files_per_subset     = DEFAULT_FILES_PER_SUBSET,

                           pythia_tmp_output    = DEFAULT_PYTHIA_OUTPUT_FILENAME,
                           output_file_index_offset = DEFAULT_OUTPUT_FILE_INDEX_OFFSET,

                           subset_index         = DEFAULT_SUBSET_INDEX,

                           upload               = False,
                           proxy                = None,
                           db_entries           = '',
                           ):
   global USE_PYTHIA_6
   
   try:
      if len(input_basenames) <= 0:
         logger.error('No inputs provided')
         raise Exception('No Input Files provided.')

      logger.info('Input basenames: \n' + '\n'.join(x for x in input_basenames))
      
      output_basename = collection_name + '._'

      logger.info('Output basenames:                  ' + str(output_basename))
      logger.info('Number of final showered events:   ' + str(nevts_shwr_final))
      logger.info('Collection Name:                   ' + str(collection_name))
      logger.info('Subset Index:                      ' + str(subset_index))
      logger.info('First Event:                       ' + str(first_event))
      logger.info('First File:                        ' + str(first_file))
      logger.info('Pythia Binary:                     ' + str(pythia_binary))
      logger.info('Pythia Inclusive:                  ' + str(pythia_inclusive))
      logger.info('Pythia Events to run:              ' + str(pythia_nevts))
      logger.info('Pythia xmldoc path:                ' + str(pythia_xmldoc_path))
      logger.info('Unweighted events per file:        ' + str(unw_evts_per_file))
      logger.info('Showered events per file:          ' + str(shwr_evts_per_file))
      logger.info('Files Per Subset:                  ' + str(files_per_subset))
      logger.info('Output File Offset:                ' + str(output_file_index_offset))
      logger.info('Upload to Grid:                    ' + str(upload))
      logger.info('Grid Proxy:                        ' + str(proxy))
      logger.info('DB Entries:                        ' + str(db_entries))

      if upload and not os.path.exists(proxy):
         raise Exception('Upload is set but proxy file unspecificed or does not exist')
      

      if pythia_binary != DEFAULT_PYTHIA_BINARY_8 and pythia_binary != DEFAULT_PYTHIA_BINARY_6:
         logger.error("Unrecognized Pythia binary: " + str(pythia_binary))
         raise Exception("Unrecognized Pythia binary: " + str(pythia_binary))

      if pythia_binary == DEFAULT_PYTHIA_BINARY_8:
         USE_PYTHIA_6 = False
      elif pythia_binary == DEFAULT_PYTHIA_BINARY_6:
         USE_PYTHIA_6 = True

      # extract text input_basenames and make a list
      for i in range(len(input_basenames)):
         input_basenames[i] = os.path.abspath(input_basenames[i])

      # record current working directory
      CWD = os.getcwd()


      # if the collection folder exists
      if os.path.exists(collection_name):
         # move to collection directory
         os.chdir(collection_name)
      else:
         # otherwise create it
         os.mkdir(collection_name)
         # move to collection directory
         os.chdir(collection_name)

      # if the subset_index is set to the default, test to see if existing subsets have already been
      # created and extract their subset index and set the subset_index and output_file_index_offset
      if subset_index == DEFAULT_SUBSET_INDEX:

         logger.info('collection folder (' + str(collection_name) + ') already exists. Now I will look for existing subsets.')
         subsets = sorted(glob.glob(collection_name+'_i*'))
         logger.info('found ' + str(len(subsets)) + ' subsets in collection folder')
         if len(subsets) > 0:
            last_subset = subsets[-1]
            index = last_subset.rfind('_i')
            extension_number = int(last_subset[index+2:index+3])
            # set the subset number to the latest extension number
            # this will be incremented by one for the next subset
            subset_index = extension_number
            logger.info('found subset number: ' + str(subset_index))
            # now find the last file index in the subset
            subset_files = sorted(glob.glob(last_subset+'/'+collection_name+'._*.tar.gz'))
            logger.info('found ' + str(len(subset_files)) + ' subset files.')
            if len(subset_files) > 0:
               last_subset_file = subset_files[-1]
               index = last_subset_file.rfind('._')
               output_file_index_offset = int(last_subset_file[index+2:index+7])
               logger.info('found subset last file index: ' + str(output_file_index_offset))
               # now look at tar ball and determine the events per file to avoid pythia
               # use the second to last file though, because the last file may not be completely full since it was 
               # just the leftovers at the end of another reorganization
               second_to_last_subset_file = collection_name+('._%05i'%(output_file_index_offset-1))
               os.system('tar zxf ' + os.path.join(last_subset,second_to_last_subset_file+'.tar.gz'))
               unwPar = AlpgenUnwParFile.AlpgenUnwParFile.read_file(second_to_last_subset_file+'.dat')
               unw_evts_per_file = unwPar.event_count
               logger.info('found unweighted events per file should be: ' + str(unw_evts_per_file))
               # remove tar files
               os.system('rm '+second_to_last_subset_file+'.{dat,events}')
            else:
               raise Exception('Did not find any subset files. Could indicate a previously failed run, so not continuing')
         else:
            # set subset to 0 for starting
            subset_index = 0
      else:
         logger.info('using user defined subset index and file offset.')
         # since I add 1 to the subset_index to handle the difference between array 
         # counting in programming that starts at 0
         # and that the indices start at 1 in rucio
         if subset_index >= 0:
            subset_index = subset_index - 1
      
      
      # get the number of events to include per file
      # this runs Pythia over a subset of events and calculates the showering efficiency
      if unw_evts_per_file == DEFAULT_UNW_EVTS_PER_FILE:
         try:
            unw_evts_per_file = get_unw_evts_per_file(
                                 input_basenames[0],
                                 shwr_evts_per_file,
                                 pythia_binary,pythia_inclusive,
                                 pythia_nevts,pythia_xmldoc_path,
                                 pythia_tmp_output)
         except Exception,e:
            logger.error('Exception running pythia: ' + str(e))
            raise
      
      logger.info(' Number of Unweighted Events per Output File: ' + str(unw_evts_per_file))
      logger.info(' Target Number of Showered events in output File: ' + str(shwr_evts_per_file))
      logger.info(' Showering Efficiency: ' + str(shwr_evts_per_file*100./unw_evts_per_file) + '%')
      logger.info(' Files per Subset: ' + str(files_per_subset))
      logger.info(' Target Showered events: ' + str(nevts_shwr_final))
      #num_subsets = nevts_shwr_final / (5000. * files_per_subset) 
      #logger.info(' Number of Subsets: ' + str(num_subsets))

      
      if files_per_subset > MAX_FILES_PER_SUBSET:
         raise MaxFilesPerSubsetExceeded(' number of files per subset exceeds the max. MAX = ' + str(MAX_FILES_PER_SUBSET) + ', user setting = ' + str(files_per_subset) + '. Try increasing the number of events per subset.')


      

      # keep track of total number of events written
      total_events = 0
      total_files = 0

      # the data block size that will be used
      read_write_data_size = DEFAULT_READ_WRITE_BLOCK_SIZE

      # keep track of all subset names
      subset_names = []
      # current subset object
      subset = None


      # loop over input basenames and reorganize
      for basename in input_basenames:
         logger.info('processing (' + str(input_basenames.index(basename)) + ' of ' + str(len(input_basenames)) + ' : ' + str(basename))

         # open input files
         input_files = InputFiles(basename)


         # create a new subset if it doesn't already exist
         if subset is None:
            subset_name = '%s_i%i1' % (collection_name,subset_index+1)
            subset_names.append(subset_name)
            subset = Subset(
                            subset_name,
                            unw_evts_per_file,
                            shwr_evts_per_file,
                            files_per_subset,
                            output_basename,
                            input_files.input_par_file,
                            pythia_inclusive,
                            output_file_index_offset,
                           )
         # calculate the data size to read
         read_write_data_size = (int(DEFAULT_READ_WRITE_BLOCK_SIZE/input_files.bytes_per_event)+1)*input_files.bytes_per_event
         # the read/write data size should not be more than the max file size or it will cause problems
         if read_write_data_size > subset.output_file_max_size:
            read_write_data_size = subset.output_file_max_size

         logger.debug('data_size = ' + str(read_write_data_size))

         data_block = input_files.read_block(read_write_data_size)
         
         while len(data_block) > 0:
            #logger.debug('while loop ' + str(len(data_block)))

            leftovers = subset.write_block(data_block)
            
            # if the subset was closed and some data left over, add it to next subset
            if leftovers is not None:
               logger.debug('in leftovers: ' + str(len(leftovers)))
               total_files += subset.output_file_counter
               subset = None
               subset_index += 1
               if num_subsets >= MAX_SUBSETS:
                  raise MaxNumSubsetExceeded(' number of subsets to produce exceeds the max. MAX = ' + str(MAX_SUBSETS) + ', user setting = ' + str(num_subsets) + '. Try increasing the number of events per subset.')

               subset_name = '%s_i1%i' % (collection_name,subset_index+1)
               subset_names.append(subset_name)
               subset = Subset(
                               subset_name,
                               unw_evts_per_file,
                               shwr_evts_per_file,
                               files_per_subset,
                               output_basename,
                               input_files.input_par_file,
                               pythia_inclusive,
                               output_file_index_offset + subset_index*files_per_subset,
                              )
               subset.write_block(leftovers)

            total_events += read_write_data_size/subset.bytes_per_event

            # check if we've written enough data
            if int(total_events/unw_evts_per_file)*5000 >= nevts_shwr_final:
               logger.debug(' total_events exceeded')
               break

            # read next data block
            data_block = input_files.read_block(read_write_data_size)


      # close the last subset
      subset.close()
      total_files += subset.output_file_counter
      subset = None
      subset_index += 1
      

      os.chdir(CWD)

      logger.info('Done creating files ' + collection_name + '  ' + str(total_events) + ' events written across ' + str(subset_index) + ' subsets with a total of ' + str(total_files) + ' files with each file containing ' + str(unw_evts_per_file) + ' events.')

      
      if upload:
         
         for subset_name in subset_names:
            rucio_bulk_upload.rucio_bulk_upload(
                     glob.glob(os.path.join(collection_name,subset_name,'*')),
                     subset_name,
                     collection_name,
                     username='usatlas2',
                     proxyfile=proxy,
                     db_entries=db_entries,
                     target_events_total=nevts_shwr_final
                     )


   except Exception,e:
      logger.exception(' received exception while running: ' + str(e))
      raise
         


class Subset:
   def __init__(self,
                subset_name,
                unw_evts_per_file,
                shwr_evts_per_file,
                files_per_subset,
                output_basename,
                input_par_file,
                pythia_inclusive,
                file_index_offset = 0,
               ):
      self.subset_name           = subset_name

      self.unw_evts_per_file     = unw_evts_per_file
      self.shwr_evts_per_file    = shwr_evts_per_file
      self.files_per_subset      = files_per_subset

      self.output_basename       = output_basename
      self.input_par_file        = input_par_file
      self.bytes_per_event       = self.input_par_file.get_bytes_per_event()

      self.output_files          = None
      self.file_index_offset     = file_index_offset
      self.output_file_counter   = 0
      self.output_file_max_size  = self.unw_evts_per_file * self.bytes_per_event
      logger.debug('output_file_max_size: ' + str(self.output_file_max_size))

      self.pythia_inclusive      = pythia_inclusive

      self.events_written        = 0

      self.closed                = False

      self.create_folder()
      pass

   def create_folder(self):
      # if no basename provided, work in current directory
      if self.subset_name is not None:
         # create folder
         os.mkdir(self.subset_name)

   def write_block(self,data_block):
      # if this subset has not been closed
      if not self.closed:
         # if no output file exists, create one
         if self.output_files is None:
            self.output_files = OutputFiles(os.path.join(self.subset_name,self.output_basename),
                                            self.output_file_counter+self.file_index_offset,
                                            self.bytes_per_event)

            logger.debug('unw_output_bytes written: ' + str(self.output_files.unw_output_bytes))

         # if data_block will overflow the output file, write small bit, close file, open new one
         if self.output_files.unw_output_bytes + len(data_block) > self.output_file_max_size:
            logger.debug(' data_block exceeds remaining file size' )
            bytes_left_to_write = self.output_file_max_size-self.output_files.unw_output_bytes
            
            self.output_files.write_block(data_block[0:bytes_left_to_write])
            self.output_files.close(self.input_par_file,self.pythia_inclusive)
            self.output_files = None
            # increment counter
            self.output_file_counter += 1
            logger.debug('output_file_counter: ' + str(self.output_file_counter))
            if self.output_file_counter >= self.files_per_subset:
               logger.debug(' output file count exceeded for subset ' )
               self.closed = True
               self.events_written += (len(data_block)-bytes_left_to_write)/self.bytes_per_event
               return data_block[bytes_left_to_write:]

            # create new file
            self.output_files = OutputFiles(os.path.join(self.subset_name,self.output_basename),
                                            self.output_file_counter+self.file_index_offset,
                                            self.bytes_per_event)
            # write remaining data to new file
            self.output_files.write_block(data_block[bytes_left_to_write:])
            self.events_written += len(data_block)/self.bytes_per_event
         else:
            self.output_files.write_block(data_block)
      else:
         logger.warning(' trying to write to a closed subset. ')
      return None

   def close(self):
      self.output_files.close(self.input_par_file,self.pythia_inclusive)
      self.output_file_counter += 1
      self.closed = True
      self.output_files = None

class OutputFiles:
   def __init__(self,output_basename,output_file_counter,bytes_per_event):
      self.output_basename          = output_basename
      self.output_file_counter      = output_file_counter
      # need a +1 everywhere because grid datasets start counting at 1!! stupid...
      self.output_unw_filename      = get_output_filename(output_basename,output_file_counter+1,'.unw')
      self.output_par_filename      = get_output_filename(output_basename,output_file_counter+1,'_unw.par')

      self.output_tarball_dat       = get_output_filename(output_basename,output_file_counter+1,'.dat')
      self.output_tarball_events    = get_output_filename(output_basename,output_file_counter+1,'.events')
      self.output_tarball           = get_output_filename(output_basename,output_file_counter+1,'.tar.gz')

      self.unw_output_bytes         = 0
      self.bytes_per_event          = bytes_per_event

      self.output_unw_file          = AlpgenUnwFile_v2.AlpgenUnwFile(self.output_unw_filename,self.bytes_per_event)

   def write_block(self,data_block):
      self.unw_output_bytes += len(data_block)
      self.output_unw_file.write(data_block)

   def close(self,par_file,pythia_inclusive):
      # update Par file settings and write it.
      new_par_file = copy.deepcopy(par_file)
      
      new_par_file.event_count = int(self.unw_output_bytes/self.bytes_per_event)
      new_par_file.lumi = new_par_file.event_count/new_par_file.cross_section

      if pythia_inclusive:
         new_par_file.add_pythia6_settings(jet_treatment=new_par_file.INCLUSIVE_JET_MATCHING)
      else:
         new_par_file.add_pythia6_settings(jet_treatment=new_par_file.EXCLUSIVE_JET_MATCHING)

      new_par_file.write_file(self.output_par_filename)

      self.output_unw_file = None

      self.create_tarball(True)

   def create_tarball(self, delete_input_files=False):
      if delete_input_files:
         shutil.move(self.output_par_filename,self.output_tarball_dat)
         shutil.move(self.output_unw_filename,self.output_tarball_events)
      else:
         shutil.copyfile(self.output_par_filename,self.output_tarball_dat)
         shutil.copyfile(self.output_unw_filename,self.output_tarball_events)
      os.system('tar zcf ' + str(self.output_tarball) + ' -C ' + os.path.dirname(self.output_tarball_dat) + ' ' + os.path.basename(self.output_tarball_dat) + ' ' + os.path.basename(self.output_tarball_events))
      if delete_input_files:
         os.remove(self.output_tarball_events)
         os.remove(self.output_tarball_dat)





class InputFiles:
   def __init__(self,basename):
      self.basename = basename

      self.input_par_filename = str(basename) + '_unw.par'
      # check that the par file exists
      if not os.path.exists(self.input_par_filename):
         logger.error('Unweighted par file does not exists: ' + self.input_par_filename)
         raise Exception('Unweighted par file does not exists')
      
      # read in unweighted parameter file
      self.input_par_file = AlpgenUnwParFile.AlpgenUnwParFile.read_file(self.input_par_filename)

      
      self.bytes_per_event = self.input_par_file.get_bytes_per_event()
      logger.info(' process: ' + self.input_par_file.get_process() + ' njets: ' + str(self.input_par_file.get_njets()) + ' bytes: ' + str(self.bytes_per_event))

      # get unweighted filename
      self.input_unw_filename = str(basename) + '.unw'
      self.unw_file_is_gzip = False
      # check that data file exists
      if not os.path.exists(self.input_unw_filename):
         if not os.path.exists(str(self.input_unw_filename) + '.gz'):
            logger.error('Unweighted event file does not exist: ' + str(self.input_unw_filename) + '(.gz)')
            raise Exception('Unweighted events file does not exists')
         else:
            self.input_unw_filename += '.gz'
            self.unw_file_is_gzip = True


      # Open unweighted events file
      self.input_unw_file = AlpgenUnwFile_v2.AlpgenUnwFile(self.input_unw_filename,self.bytes_per_event)

   def read_block(self,block_byte_size):
      if self.input_unw_file.file is not None:
         return self.input_unw_file.file.read(block_byte_size) # returns < block_byte_size when EOF reached


def get_output_filename(basename,counter,extension):
   return '%s%05i%s' %(basename,counter,extension)


def get_unw_evts_per_file(
                           input_basename,
                           showered_evts_per_file,
                           binary = DEFAULT_PYTHIA_BINARY,
                           inclusive = DEFAULT_PYTHIA_JET_MATCHING_INCLUSIVE,
                           shower_nevts = DEFAULT_PYTHIA_NEVTS,
                           xmldoc_path = DEFAULT_PYTHIA_XMLDOC,
                           output_filename = DEFAULT_PYTHIA_OUTPUT_FILENAME,
                         ):
   PYTHIA_SEARCH_STRING_8 = '| sum                                                |'
   PYTHIA_SEARCH_STRING_6 = 'I   0 All included subprocesses    I'
   if USE_PYTHIA_6:
      PYTHIA_SEARCH_STRING = PYTHIA_SEARCH_STRING_6
   else:
      PYTHIA_SEARCH_STRING = PYTHIA_SEARCH_STRING_8
   
   cmd = binary
   if inclusive:
      cmd += ' -i' # turn on inclusive jets
   cmd += ' -n ' + str(shower_nevts) # set number of events to process
   if not USE_PYTHIA_6:
      cmd += ' -x ' + str(xmldoc_path) # set xmldoc path
   cmd += ' -o ' + str(output_filename) # set output filename (to be deleted later)
   cmd += ' -a ' + str(input_basename) # add input file basename
   logger.info("running pythia: " + cmd)
   p = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   pythia_tried_events = 0
   pythia_selected_events = 0
   pythia_accepted_events = 0
   pythia_efficiency = 0.
   for line in p.stdout:
      logger.info(line[0:-1])
      if PYTHIA_SEARCH_STRING in line:
         split = line.split()
         if USE_PYTHIA_6:
            pythia_tried_events = int(split[7])
            pythia_accepted_events = int(split[6])
            pythia_selected_events = pythia_accepted_events
         else:
            pythia_tried_events = int(split[3])
            pythia_selected_events = int(split[4])
            pythia_accepted_events = int(split[5])
         if pythia_accepted_events != shower_nevts:
            logger.warning(' Number of pythia events accepted = ' + str(pythia_accepted_events) + ' when ' + str(shower_nevts) + ' were requested. Pythia tried ' + str(pythia_tried_events) + '.')
         pythia_efficiency = float(pythia_accepted_events)/float(pythia_tried_events)
         logger.info(' Pythia accepted ' + str(pythia_accepted_events) + ' events out of ' + str(pythia_tried_events) + ' for efficiency of ' + str(pythia_efficiency))
         


   p.wait()
   if p.returncode != 0:
      logger.error('pythia exited with non-zero return code: ' + str(p.returncode) + '. STDERR = ')
      for line in p.stderr:
         logger.error(line[0:-1])
      if os.path.exists(output_filename):
         os.remove(output_filename)
      logger.error('Pythia Exited with non-zero return code')
      raise Exception('Pythia Exited with none-zero return code')

   # remove pythia file
   os.remove(output_filename)

   # calculate the number of unw events per file rounded to the nearest 1000
   unw_evt_per_file = int( (float(showered_evts_per_file) / float(pythia_efficiency) / 1000.) + 1 ) * 1000
   logger.info(' unweighted events per file:  %i' % unw_evt_per_file)
   logger.info(' pythia efficiency: %10.8e' % pythia_efficiency)

   return unw_evt_per_file


if __name__ == "__main__":
   sys.exit(main())
