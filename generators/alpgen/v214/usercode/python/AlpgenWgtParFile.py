from AlpgenCmdFile import AlpgenCmdFile
from FileTools import GetLineCount
import subprocess,os,shutil,logging
logger = logging.getLogger(__name__)




def replace_event_count(input_filename):
   input_file = AlpgenCmdFile()
   input_file.read(input_filename)
   filename_base = input_file.filename_base

   par_filename = filename_base + '.par'
   weighted_events_filename = filename_base + '.wgt'
   new_par_filename = filename_base + '.par.new'

   # read the number of events in the weighted event file
   new_event_count = GetLineCount(weighted_events_filename)

   update_event_count(par_filename,new_par_filename,new_event_count)

def update_event_count(input_filename,output_filename,new_event_count):

   par_file = open(input_filename,'r')
   new_par_file = open(output_filename,'w')

   # find the line where the number of weighted events are located
   numbers_are_next = False
   for line in par_file:
      if line.find('number wgted evts in the file') >= 0:
         # write current line to new file
         numbers_are_next = True
      elif numbers_are_next:
         # get next line with numbers on it, format is "# wgt events"  "sigma" "error" "maxwgt"
         number_strings = line.split()
         sigma = number_strings[1]
         if sigma[0] is '.':
            sigma = '0' + sigma
         error = number_strings[2]
         if error[0] is '.':
            error = '0' + error
         maxwgt = number_strings[3]
         if maxwgt[0] is '.':
            maxwgt = '0' + maxwgt
         new_line = ("%15.1f %s %s %s\n" % (new_event_count,sigma,error,maxwgt))
         logger.debug('old par line: ' + line)
         logger.debug('new par line: ' + new_line)
         new_par_file.write(new_line)
         numbers_are_next = False
         continue
      new_par_file.write(line)

   par_file.close()
   new_par_file.close()

