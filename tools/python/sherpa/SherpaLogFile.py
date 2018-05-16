#!/usr/bin/env python
import os,sys,optparse,logging,datetime
logger = logging.getLogger(__name__)

class SherpaLogFile:
   RUN_CMD_STR = 'Running command:'
   START_TIME_STR = 'The local time is'
   FIRST_EVENT_STR = 'Event 1 ('
   FINAL_TIME_STR = 'Time:'
   NUMBER_OF_THREADS_STR = 'Process_Group::CalculateTotalXSec(): Set number of threads'
   NUMBER_OF_RANKS_STR = 'My_MPI::SetUpSendRecv(): Analyzing MPI environment {'

   def __init__(self,filename):
      self.filename = filename
      self.start_time = None
      self.command = None
      self.final_time = None
      self.evt_gen_sec = None
      self.init_sec = None
      self.num_threads = 1
      self.num_ranks   = 1
      self.num_events_per_rank = None

   def __str__(self):
      s = '\n'
      s += 'Contents of Sherpa Log File: ' + str(self.filename) + '\n'
      if self.command: s += ' command: ' + self.command + '\n'
      if self.start_time: s += ' start time: ' + self.start_time.strftime('%c') + '\n'
      if self.final_time: s += ' final time: ' + self.final_time.strftime('%c') + '\n'
      if self.evt_gen_sec: s += ' time in event generation: ' + str(self.evt_gen_sec) + ' seconds\n'
      if self.init_sec: s += ' time in initialize: ' + str(self.init_sec) + ' seconds \n'
      s += ' number of threads: ' + str(self.num_threads) + '\n'
      s += ' number of ranks: ' + str(self.num_ranks) + '\n'
      if self.num_events: s+= 'number of events: ' + str(self.num_events) + '\n'
      s += '\n'
      return s

   @staticmethod
   def read_file(filename):
      log = SherpaLogFile(filename)
      try:
         file = open(filename)
      except:
         logger.exception('Exception received while opening file ' + filename)
         raise

      in_ranks = False
      for line in file:
         #line = line[0:-1]
         line = line.replace('\r','')
         line = line.replace('\n','')
         #logger.info(line)
         if   SherpaLogFile.RUN_CMD_STR in line:
            log.command = line[len(SherpaLogFile.RUN_CMD_STR)+1:].strip()
         elif SherpaLogFile.START_TIME_STR in line:
            log.start_time = datetime.datetime.strptime(line[len(SherpaLogFile.START_TIME_STR) + 1:-1],'%c')
         #elif SherpaLogFile.FIRST_EVENT_STR in line:
            #log.first_event_time = datetime.datetime.strptime(line[line.find('ETA:')])
         elif SherpaLogFile.FINAL_TIME_STR in line:
            log.final_time = datetime.datetime.strptime(line[line.find('on ')+3:],'%c')
         elif 'Event' in line and 's total )' in line:
            # extract run time
            end_index = line.find('s total )')
            begin_index = line.rfind('(',0,end_index)
            log.evt_gen_sec = int(line[begin_index + 1:end_index-1])
            # extract total events
            end_index = begin_index
            begin_index = line.rfind('Event',0,end_index) + len('Event')
            log.num_events_per_rank = int(line[begin_index+1:end_index-1])
         elif SherpaLogFile.NUMBER_OF_THREADS_STR in line:
            log.num_threads = int(line[len(SherpaLogFile.NUMBER_OF_THREADS_STR):].replace('.','').strip())
         elif SherpaLogFile.NUMBER_OF_RANKS_STR in line:
            in_ranks = True
         elif in_ranks:
            if 'Rank' in line:
               log.num_ranks += 1
            else:
               in_ranks = False
               log.num_ranks -= 1





      total_time = log.final_time - log.start_time
      log.init_sec = total_time.seconds - log.evt_gen_sec

      return log


def main():
   logging.basicConfig(level=logging.INFO)

   parser = optparse.OptionParser(description='Test SherpaLogFile class')
   parser.add_option('-i','--input',dest='input',help='input')
   options,args = parser.parse_args()

   if options.input is None:
      parser.error('Must specify -i')

   log = SherpaLogFile.read_file(options.input)

   logger.info(str(log))


   


if __name__ == "__main__":
   main()
