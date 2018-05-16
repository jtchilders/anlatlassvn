#!/usr/bin/env python

import os,sys,subprocess,optparse,datetime,glob

CONDOR_LOG='condor_log.txt'
COBALT_LOG_GLOB='*.cobaltlog'
ALPGEN_UNW='alpout.unw'
ALPGEN_UNW_PAR='alpout_unw.par'
ALPGEN_POSTSUBMIT='alpgen_postsubmit.sh'
ALPGEN_INPUT_IMODE1 ='alpout.input.1'


def main():
   parser = optparse.OptionParser(description='extract spreadsheet information from a job')
   parser.add_option('-i','--input-folder',dest='inputFolder',help='Folder holding the job files')
   options,args = parser.parse_args()

   if options.inputFolder is None:
      parser.error('Must specify an input folder using -i')

   
   print ' getting information from folder: ' + str(options.inputFolder)
   
   get_input_file_info(options.inputFolder)

   get_condor_runtime(options.inputFolder)

   get_cobalt_queue_run_time(options.inputFolder)

   get_num_unw_evt(options.inputFolder)

   get_unw_file_size(options.inputFolder)


def get_condor_runtime(folder):
   filename = os.path.join(folder,CONDOR_LOG)
   try:
      file = open(filename,'r')
   except:
      print 'Error opening file: ' + str(filename)
      return
   
   for line in file:
      if line.find('Run Remote Usage') >= 0:
         words = line.split()
         time = words[2][:-1]
         print '%50s %25s %25s' % (filename,' Condor Run Time: ',time)
         return
   
   
def get_cobalt_queue_run_time(folder):
   # take the most recent cobaltlog
   filename = glob.glob(folder + '/' + COBALT_LOG_GLOB)[-1]
   try:
      file = open(filename,'r')
   except:
      print 'Error opening file: ' + str(filename)
      return None
   
   queue_start_time = None
   queue_end_time   = None
   run_start_time   = None
   run_end_time     = None
   for line in file:
      #print line[:-1]
      time = line[0:24]
      timezone = line[25:31]
      if line.find('submitted with cwd set to:') >= 0:
         # Fri Jul 25 21:09:42 2014 -0500 (CDT)
         queue_start_time = datetime.datetime.strptime(time,'%c') - get_timezone_offset(timezone)
         continue
      elif queue_start_time is not None and queue_end_time is None:
         queue_end_time = datetime.datetime.strptime(time,'%c') - get_timezone_offset(timezone)
         run_start_time = queue_end_time
         continue
      elif line.find('initiating job cleanup and removal') >= 0:
         run_end_time = datetime.datetime.strptime(time,'%c') - get_timezone_offset(timezone)
         continue
   
   if queue_start_time is not None and queue_end_time is not None and run_end_time is not None:
      print '%65s' % filename
      print '     %25s %25s %25s %25s ' % ('queue start','queue end','run start','run end')
      print '     %25s %25s %25s %25s ' % (str(queue_start_time),str(queue_end_time),str(run_start_time),str(run_end_time))
      print '     %25s %25s %25s %25s ' % ('queue time:',str(queue_end_time-queue_start_time),
                                           'run time:',str(run_end_time-run_start_time))


def get_timezone_offset(timezone):
   # time zone is -0500 for central or -0400 for eastern, etc.
   timezone = int(timezone)
   timezone = int(timezone / 100)
   return datetime.timedelta(hours=timezone)

def get_num_unw_evt(folder):
   filename = os.path.join(folder,ALPGEN_UNW_PAR)
   try:
      file = open(filename,'r')
   except:
      print 'Error opening file: ' + str(filename)
      return None
   
   for line in file:
      if line.find('! unwtd events, lum (pb-1)') >= 0:
         str_num_evts = line.split()[0]
         num_evts = int(str_num_evts)
         print '%65s %25s %25d' % (filename,'number unw events',num_evts)
         return

def get_unw_file_size(folder):
   filename = os.path.join(folder,ALPGEN_UNW)
   p = subprocess.Popen(['ls','-l',filename],stderr=subprocess.PIPE,stdout=subprocess.PIPE)
   stdout,stderr = p.communicate()
   if p.returncode != 0:
      print ' ERROR getting unw file size'
      return
   size = int(stdout.split()[4])
   print '%65s %25s %25d' % (filename,'unw file size',size)

def get_input_file_info(folder):
   filename = os.path.join(folder,ALPGEN_INPUT_IMODE1)
   
   iter_events = 0
   iters = 0
   events = 0
   iseed1 = 0
   iseed2 = 0
   iseed3 = 0
   iseed4 = 0
   for line in open(filename):
      if line.find('Num events/iteration') >= 0:
         x = line.split()
         iter_events = int(x[0])
         iters = int(x[1])
      elif line.find('Num events generated') >= 0:
         x = line.split()
         events = int(x[0])
      elif line.find('iseed1') >= 0:
         iseed1 = int(line.split()[1])
      elif line.find('iseed2') >= 0:
         iseed2 = int(line.split()[1])
      elif line.find('iseed3') >= 0:
         iseed3 = int(line.split()[1])
      elif line.find('iseed4') >= 0:
         iseed4 = int(line.split()[1])

   print 'evts/iter = %10i/%3i evts = %10i' % (iter_events,iters,events)
   print 'iseed1/2/3/4 = %15i %15i %15i %15i' % (iseed1,iseed2,iseed3,iseed4)






   

if __name__ == '__main__':
   main()


