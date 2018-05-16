#!/usr/bin/env python
import os,sys,optparse,logging,glob,json,decimal
import AlpgenInputFile
from mysql import mysql_wrapper
import CalcMean
logger = logging.getLogger(__name__)

def main():
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='')
   parser.add_option('-g','--glob',dest='glob',help='glob')
   parser.add_option('-e','--do-evtgen',dest='do_evtgen',action='store_true',default=False,help='do the event generation analysis.')
   parser.add_option('-f','--do-gridgen',dest='do_gridgen',action='store_true',default=False,help='do the grid generation analysis.')
   options,args = parser.parse_args()

   
   manditory_args = [
                     'glob',
                  ]

   for man in manditory_args:
      if options.__dict__[man] is None:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)
   
   #mysql = mysql_wrapper.MySQL(db_server_password='1w@ANLHEP')

   file_list = glob.glob(options.glob)

   unique_file_list = {}

   for filename in file_list:
      #logger.info('filename: ' + filename)
      basename = os.path.basename(filename)
      if basename in unique_file_list:
         pass
      else:
         unique_file_list[basename] = filename

   if options.do_gridgen:
      analyze_gridgen(unique_file_list)
   elif options.do_evtgen:
      analyze_evtgen(unique_file_list)
   else:
      parser.print_help()
      sys.exit(-1)



def analyze_gridgen(unique_file_list):

   rates = {}
   for basename,filename in unique_file_list.iteritems():
      #logger.info('filename: %100s' % filename)

      imode0 = False
      gridmode0 = False
      n_evts_per_itr = 0
      n_itr = 0
      n_evts = 0
      run_time = 0
      inline_process = ''
      with open(filename) as fp:
         line = fp.readline()
         while line != '':
            #logger.info(line[:-1])
            if 'imode =' in line:
               parts = line.split()
               if int(parts[2]) == 0:
                  imode0 = True
            elif imode0 and 'grid source =' in line:
               parts = line.split()
               if int(parts[3]) == 0:
                  gridmode0 = True
            elif imode0 and gridmode0 and 'N(events)/iterations and N(iterations) =' in line:
               parts = line.split()
               n_evts_per_itr = int(float(parts[4]))
               n_itr = int(float(parts[6]))
            elif imode0 and gridmode0 and 'number of events to generate =' in line:
               parts = line.split()
               n_evts = int(float(parts[6]))
            elif imode0 and gridmode0 and 'FINALIZE RANK 00000000 AT' in line:
               parts = line.split()
               run_time = int(float(parts[4]))
            elif 'kperp is then rescaled by ktfac' in line:
               line = fp.readline()
               while line.strip() == '':
                  line = fp.readline()
               if len(line) > 0:
                  inline_process = line[:-1]
                  #logger.info('process = ' + process + ' ' + line)

            #input('here')
            line = fp.readline()
      process = get_process(inline_process) + ' + ' + str(get_njets(inline_process)) + ' jets'

      # calculate the rate if I have all the values I need
      #logger.info(str([n_evts_per_itr,n_itr,n_evts,run_time,len(process)]))
      if all(x != 0 for x in [n_evts_per_itr,n_itr,n_evts,run_time,len(process)]):
         rate = float(n_evts_per_itr*n_itr + n_evts)/float(run_time)

         # store the rate by process name
         if process in rates:
            rates[process].add_value(rate)
         else:
            m = CalcMean.CalcMean()
            m.add_value(rate)
            rates[process] = m
         #logger.info('    %15i %5i %15i %15i  %s' % (n_evts_per_itr,n_itr,n_evts,run_time,process))

   
   logger.info('\n' + json.dumps(rates, sort_keys=True,indent=3, separators=(',', ': '),cls=CalcMean.CalcMeanEncoder))
   #logger.info('\n' + json.dumps(rB, sort_keys=True,indent=3, separators=(',', ': ')))


def analyze_evtgen(unique_file_list):

   rates = {}
   for basename,filename in unique_file_list.iteritems():
      logger.debug('filename: %100s' % filename)

      # analyze the Alpgen Log File
      imode1 = False
      gridmode2 = False
      n_evts_per_itr = 0
      n_itr = 0
      n_evts = 0
      run_time = 0
      process = ''
      failed = False
      evgen_start_sec = 0
      unw_start_sec = 0
      agg_start_sec = 0
      job_end_sec = 0
      try:
         with open(filename) as fp:
            line = fp.readline()
            while line != '':

               if 'imode =' in line:
                  parts = line.split()
                  if int(parts[2]) == 1:
                     imode1 = True
               elif imode1 and 'grid source =' in line:
                  parts = line.split()
                  if int(parts[3]) == 2:
                     gridmode2 = True
               elif imode1 and gridmode2 and 'N(events)/iterations and N(iterations) =' in line:
                  parts = line.split()
                  n_evts_per_itr = int(float(parts[4]))
                  n_itr = int(float(parts[6]))
               elif imode1 and gridmode2 and 'number of events to generate =' in line:
                  parts = line.split()
                  n_evts = int(float(parts[6]))
               elif imode1 and gridmode2 and 'FINALIZE RANK 00000000 AT' in line:
                  parts = line.split()
                  run_time = int(float(parts[4]))
               elif 'kperp is then rescaled by ktfac' in line:
                  line = fp.readline()
                  while line.strip() == '':
                     line = fp.readline()
                  if len(line) > 0:
                     inline_process = line[:-1]
               elif 'Failed with exit code:' in line:
                  failed = True
               elif 'Generating weighted events' in line:
                  parts = line.split()
                  evgen_start_sec = int(parts[3])
               elif 'Unweighting events' in line:
                  parts = line.split()
                  unw_start_sec = int(parts[2])
               elif 'Aggregating results' in line:
                  parts = line.split()
                  agg_start_sec = int(parts[2])
               elif 'Unweighting complete, job complete' in line:
                  parts = line.split()
                  job_end_sec = int(parts[4])

               line = fp.readline()
      except: continue

      if failed: continue


      process = get_process(inline_process) + ' + ' + str(get_njets(inline_process)) + ' jets'
      
      # analyze the paramater file to get total unw events
      '''
      par_file = os.path.join(os.path.dirname(filename),'alpout_unw.par')
      num_unw_evts = 0
      if os.path.exists(par_file):
         for line in open(par_file):
            if '! unwtd events, lum (pb-1)' in line:
               parts = line.split()
               num_unw_evts = int(parts[0])
               break
      else:
         logger.error('no par file: ' + par_file)
      '''

      # analyze the cobaltlog, get number of nodes and ranks per node
      log_file = filename.replace('output','cobaltlog')
      nodes = 0
      ranks_per_node = 0
      if os.path.exists(log_file):
         for line in open(log_file):
            if 'qsub' in line:
               parts = line.split()
               for i in range(len(parts)):
                  if parts[i] == '-n':
                     nodes = int(parts[i+1])
               ranks_per_node = int(parts[-1])
               break
      else:
         logger.error('no cobaltlog file: ' + log_file)

      
      # calculate all the rates
      if all(x != 0 for x in [n_evts,run_time,len(process),evgen_start_sec,unw_start_sec,agg_start_sec,job_end_sec]):
         # the rate from start of event gen to end of unweighting
         rate = float(n_evts)/float(job_end_sec - evgen_start_sec)

         if process in rates:
            rates[process].add_value(rate)
         else:
            m = CalcMean.CalcMean()
            m.add_value(rate)
            rates[process] = m

   logger.info('\n' + json.dumps(rates, sort_keys=True,indent=3, separators=(',', ': '),cls=CalcMean.CalcMeanEncoder))
   
   


def get_njets(process):
   if 'production' in process:
      return 0

   if '+' in process:
      parts = process.split('+')
      count = 0
      if len(parts) == 3:
         count += int(parts[2].replace('jets',''))
      else:
         count += int(parts[1].replace('jets',''))
      return count

   return 0

def get_process(process):
   if 'production' in process:
      parts = process.split()
      if 'Z' == parts[0].strip():
         return 'zjet'
      elif 'W' == parts[0].strip():
         return 'wjet'
      return ''

   elif '+' in process:
      parts = process.split('+')
      if len(parts) == 2:
         if 'Z' == parts[0].strip():
            return 'zjet'
         elif 'W' == parts[0].strip():
            return 'wjet'
         elif 'W b bbar' in parts[0] or 'W c cbar' in parts[0]:
            return 'wqq'
      elif 'W + c' in process:
         return 'wcjet'
   elif 'W b bbar' == process.strip():
      return 'wqq'
   elif 'W c cbar' == process.strip():
      return 'wqq'
   logger.info('did not identify process: ' + process)
   return ''


def get_object_count(process):
   #logger.debug(process)
   if 'production' in process:
      return 1
   elif '+' in process:
      parts = process.split('+')
      count = 0
      if len(parts) == 3:
         count += len(parts[0].split())
         count += len(parts[1].split())
         count += int(parts[2].replace('jets',''))
      else:
         count += len(parts[0].split())
         count += int(parts[1].replace('jets',''))
      return count
   else:
      parts = process.strip().split()
      return len(parts)


if __name__ == "__main__":
   main()
