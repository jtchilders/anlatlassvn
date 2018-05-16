#!/usr/bin/env python
import os,sys,subprocess,optparse,logging,traceback
logger = logging.getLogger(__name__)
from AlpgenUnwParFile import AlpgenUnwParFile
from AlpgenUnwFile import AlpgenUnwFile,ALPGEN_ATHENA_FILEBASE

TMP_FOLDER='tmp_alpgen_unw_split'
TARBALL_FOLDER='tarballs'
DEFAULT_INPUTBASE='alpout'
DEFAULT_OUTPUTBASE='split'
DEFAULT_MAXFILES=-1
DEFAULT_NEVT=5000
DEFAULT_OUTPUT_OFFSET = 0
DEFAULT_FIRST_EVENT = 0
DEFAULT_FIRST_FILE = 0
DEFAULT_SKIP_LAST_FILE = False
DEFAULT_MAKE_TARBALLS = True
DEFAULT_TARBALL_ONLY = False

def main():
   logging.basicConfig(format='%(asctime)s %(levelname)s:%(name)s:%(message)s',level=logging.INFO)
   parser = optparse.OptionParser(description='Divide alpgen unweighted events file and parameters file, create tarball for each subset.')
   parser.add_option('-n','--evnt-per-file',dest='nevt',help='Number of events per file [DEFAULT='+str(DEFAULT_NEVT)+']',type='int',default=DEFAULT_NEVT)
   parser.add_option('-i','--input-base',dest='input_base',help='base filename for the input alpgen files, can be a comma separated list.[DEFAULT='+DEFAULT_INPUTBASE+'].',default=DEFAULT_INPUTBASE)
   parser.add_option('-o','--output-base',dest='output_base',help='Output filename base [DEFAULT='+DEFAULT_OUTPUTBASE+'].',default=DEFAULT_OUTPUTBASE)
   parser.add_option('-m','--max-files',dest='max_files',help='Maximum number of files to output [DEFAULT='+str(DEFAULT_MAXFILES)+'].',type='int',default=DEFAULT_MAXFILES)
   parser.add_option('-f','--output-offset',dest='output_offset',help='The offset to start with when creating numbered outputfiles. [DEFAULT='+str(DEFAULT_OUTPUT_OFFSET)+']',type='int',default=DEFAULT_OUTPUT_OFFSET)
   parser.add_option('-s','--skip-tarball',dest='makeTarballs',help='Skip the creation of gzipped tarballs for athena input.',action='store_false',default=DEFAULT_MAKE_TARBALLS)
   parser.add_option('-l','--skip-last-file',dest='skip_last_file',help='The last file created from each base filename will not reach the "max-files" event number. This flag will cause this last file to be omitted to avoid a file with a different number of events in it.',action='store_true',default=DEFAULT_SKIP_LAST_FILE)
   parser.add_option('-t','--tarball-base',dest='tarball_base',help='The basename for the tarball files.',default=ALPGEN_ATHENA_FILEBASE)
   parser.add_option('-a','--first-event',dest='first_event',help='This is essentially at what event number to begin.[DEFAULT='+str(DEFAULT_FIRST_EVENT)+']',default=DEFAULT_FIRST_EVENT,type='int')
   parser.add_option('-b','--first-file',dest='first_file',help='This sets the number of files to skip before beginning. File numbers will start with zero still. [DEFAULT='+str(DEFAULT_FIRST_FILE)+']',default=DEFAULT_FIRST_EVENT,type='int')
   parser.add_option('-c','--output-path',dest='output_path',help='The path where the split files are written. [default = "' + TMP_FOLDER + '"]',default=TMP_FOLDER)
   parser.add_option('--tarball-only',dest='tarball_only',help='Usually the tarballs are in a subdirectory. For datasets you might want the tarballs only and not the .unw and _unw.par files. Turn this flag on to get just the tarballs.',default=DEFAULT_TARBALL_ONLY,action='store_true')
   parser.add_option('-u','--update-par',dest='update_par',help='Enabling this flag adds the Pythia6 specific settings to the _unw.par files',action='store_true',default=False)
   parser.add_option('--jet-match-excl',dest='jet_match_excl',help='Used with "update-par" option. Sets Pythia to use exclusive jet matching',action='store_true',default=False)
   parser.add_option('--jet-match-incl',dest='jet_match_incl',help='Used with "update-par" option. Sets Pythia to use inclusive jet matching',action='store_true',default=False)
   options,args = parser.parse_args()

   jet_matching = None
   if options.update_par:
      if not options.jet_match_incl and not options.jet_match_excl:
         parser.error('must set --jet-match-excl or --jet-match-incl when using --update-par')
      elif options.jet_match_incl:
         jet_matching = AlpgenUnwParFile.INCLUSIVE_JET_MATCHING
      elif options.jet_match_excl:
         jet_matching = AlpgenUnwParFile.EXCLUSIVE_JET_MATCHING
   
   return split_alpgen_unw(
                     options.input_base.split(','),
                     options.nevt,
                     options.output_base,
                     options.max_files,
                     options.output_offset,
                     options.makeTarballs,
                     options.skip_last_file,
                     options.tarball_base,
                     options.first_event,
                     options.first_file,
                     options.output_path,
                     options.tarball_only,
                     options.update_par,
                     jet_matching,
                    )
   
def split_alpgen_unw(input_basenames,
                     nevts,
                     output_base,
                     max_files=DEFAULT_MAXFILES,
                     output_offset = DEFAULT_OUTPUT_OFFSET,
                     makeTarballs = DEFAULT_MAKE_TARBALLS,
                     skip_last_file = DEFAULT_SKIP_LAST_FILE,
                     tarball_base= ALPGEN_ATHENA_FILEBASE,
                     first_event = DEFAULT_FIRST_EVENT,
                     first_file = DEFAULT_FIRST_FILE,
                     output_path = TMP_FOLDER,
                     tarball_only = DEFAULT_TARBALL_ONLY,
                     update_par = False,
                     jet_matching = None,
                    ):
   
   # make a temporary directory to hold all the files
   try:
      if not os.path.exists(output_path):
         os.mkdir(output_path)
   except:
      logger.exception('Error trying to create folder ' + output_path + ' Exception: ' + str(sys.exc_info()[1]) + '.')
      return -1


   # divide the alpgen unweighted events/parameter file into smaller bits
   split_basenames = []
   for input_basename in input_basenames:
      num_files = len(split_basenames)
      # if all the files requested have been produced then end execution
      if max_files - num_files < 0:
         break;
      offset = num_files + output_offset
      

      alpgenUnwParFile = AlpgenUnwParFile.read_file(input_basename+'_unw.par')

      # add pythia 6 settings if requested
      if update_par:
         if jet_matching is None:
            logger.error('Must specify jet_matching when using update_par')
            raise Exception('Must specify jet_matching when using update_par')
         alpgenUnwParFile.add_pythia6_settings(jet_treatment=jet_matching)

      alpgenUnwFile = AlpgenUnwFile(input_basename+'.unw',alpgenUnwParFile)
      logger.info(' file contains ' + str(alpgenUnwParFile.event_count) + ' events')

      # if the first event number is larger than the number of events in this file, move on
      if first_event > alpgenUnwParFile.event_count:
         # reduce first_event by the number in this file so it is properly compared to the next file
         first_event = first_event - alpgenUnwParFile.event_count
         # continue to next file
         continue
      # if the first file is larger than the number of output files that can be made 
      # with this input file, then move on
      if first_file > int(alpgenUnwParFile.event_count/nevts):
         # reduce first_file by the number of files that would have been created so
         # the next file is compared properly
         first_file = first_file - int(alpgenUnwParFile.event_count/nevts)
         # continue to next file
         continue
      
      logger.info(' splitting files with basename: ' + str(input_basename))
      logger.info(' events per file: ' + str(nevts))
      logger.info(' output filename base: ' + str(output_base))
      logger.info(' max files: ' + str(max_files - num_files))
      logger.info(' output file number offset: ' + str(output_offset))
      logger.info(' first event: ' + str(first_event))
      logger.info(' first file: ' + str(first_file))
      logger.info(' tarball base name: ' +str(tarball_base))
      if makeTarballs:
         logger.info(' going to make tarballs')
      if skip_last_file:
         logger.info(' going to skip the events at the end of input files which would lead to partially filled files.')

      # split this file
      try:
         basenames =  alpgenUnwFile.SplitFile(nevts,
                                                    os.path.join(output_path,output_base),
                                                    max_files - num_files,
                                                    offset,
                                                    skip_last_file,
                                                    first_event,
                                                    first_file,
                                                   )
         split_basenames += basenames
         first_event = 0
         first_file = 0
      except:
         logger.exception('Error while splitting alpgen file ' + alpgenUnwFile.filename + ' Exception: ' + str(sys.exc_info()[1]) )
         return -2

   if makeTarballs:
      # now make the tarball of each unw/unw.par file set
      logger.info(' looping over all files and making tarball of each .unw and _unw.par set. ')
      try:
         path = os.path.join(output_path,TARBALL_FOLDER)
         if not os.path.exists(path):
            os.mkdir(path)
      except:
         logger.exception('Error trying to create folder ' + os.path.join(output_path,TARBALL_FOLDER) + ' Exception: ' + str(sys.exc_info()[1])+'.')
         return -3
      try:
         alpgenUnwFile.MakeDataSetTarball(
               split_basenames,
               os.path.join(output_path,TARBALL_FOLDER,tarball_base),
               output_offset,
             )
      except:
         logger.exception('Error trying to create tarballs Exception: ')
         traceback.print_exc()
         return -4
   
   # replace the .unw and _unw.par files with the tarballs
   if tarball_only:
      os.system('rm ' + output_path + '/*[.unw,_unw.par]')
      os.system('mv ' + os.path.join(output_path,TARBALL_FOLDER) + '/* ' + output_path + '/')
      os.system('rmdir ' + os.path.join(output_path,TARBALL_FOLDER))
   

   return len(split_basenames) 
   

if __name__ == '__main__':
   ret = main()
   if ret > 0: sys.exit(0)
   else: sys.exit(ret)


