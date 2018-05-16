#! /usr/bin/python

#####################################################################
# This replaces the original Alpgen.sh and Alpmaker.c               #
#  with a single piece of code that does everything                 #
#  necessary to run an Alpgen 4-vector job on an HPC.               #
#     input file: input                                             #
#     output files: input.0       (warmup)                          #
#                   input.1       (event generation)                #
#                   input.2       (unweighting)                     #
#                   alpgen.dag    (Condor DAG submit file)          #
#                   warmup.cmd    (Condor submit for warmup phase)  #
#                   hpcwait.cmd   (Condor submit file for HPC phase #
#                                - or fake HPC phase)               #
#                   unweight.cmd  (Condor submit file for           #
#                                  unweighting phase)               #
#                   alpgen.xml    (XML file to send to HPC          #
#                                  scheduler)                       #
#                   pdfname.txt   (Name of the PDF file;            # 
#                                  makes future scripts simpler     #
#                                                                   #
#                                                                   #
#                                                                   #
#     Returns: 0 (normal completion)                                #
#              1 (problem with parameters)                          #
#              2 (problem with input file)                          #
#              3 (problem with output files)                        #
#              200 (internal error)                                 #
#####################################################################




import math, os, sys, time, shutil
sys.path.append("/users/hpcusers/svn/jobpayloads")

from Computer import Computer,computers
from CondorJob import CondorJob



class Initialize:
# Collects options from environment and command-line
   def __init__(self):
      import getopt
   
      # Parse the command line  

      self.verbose = False
      self.testmode = False
      self.doPythia = True
      self.submit = True

      warmup_1_string = '10000000'
      warmup_2_string = '20000000'
      warmup_rounds_string = '8'
     
      self.process_name = ''
      event_number_string = ''
     
      self.machine_list = []

      try:
         options, remainder = getopt.gnu_getopt(sys.argv[1:], 
                                                '1:2:m:n:p:r:tvhsj', 
                                                ['warmup1=', 'warmup2=','machine=',
                                                 'warmuprounds=','events=','nosubmit',
                                                 'verbose','skipPythia','help',
                                                ]
                                               )

      except getopt.GetoptError:
         print 'usage: python ' + sys.argv[0] + ' [options]'
         print '              -1 or --warmup1           # evts for first warmup generation'
         print '              -r or --warmuprounds      # of rounds to run "warmup1" evts'
         print '              -2 or --warmup2           # evts for last warmup generation'
         print '              -n or --events            # evts for weighted generation'
         print '              -p or --process           alpgen process to use for generation'
         print '              -m or --machine           which machine to run on'
         print '              -j or --no-submit         write all output files, but do not submit anything. For testing'
         print '              -s or --skip-pythia       skip running pythia'
         print '              -v or --verbose           lots of text'
         print '              -h or --help              print this again'
         print '  Enjoy!'
         sys.exit(1)	  

      for opt, arg in options:
         if opt in ('-1', '--warmup1'):
             warmup_1_string = arg
         elif opt in ('-2', '--warmup2'):
             warmup_2_string = arg  
         elif opt in ('-m', '--machine'):
             self.machine_list.append(arg)
         elif opt in ('-n', '--events'):
             event_number_string = arg
         elif opt in ('-p', '--process'):
             self.process_name = arg
         elif opt in ('-r', '--warmuprounds'):
             warmup_rounds_string = arg
         elif opt in ('-s', '--skip-pythia'):
             self.doPythia = False
         elif opt in ('-j','--no-submit'):
             self.submit = False
         elif opt in ('-v', '--verbose'):
             self.verbose = True
         elif opt in ('-h','--help'):
            print 'usage: python ' + sys.argv[0] + ' [options]'
            print '              -1 or --warmup1           # evts for first warmup generation'
            print '              -r or --warmuprounds      # of rounds to run "warmup1" evts'
            print '              -2 or --warmup2           # evts for last warmup generation'
            print '              -n or --events            # evts for weighted generation'
            print '              -p or --process           alpgen process to use for generation'
            print '              -m or --machine           which machine to run on'
            print '              -v or --verbose           lots of text'
            print '              -j or --no-submit         write all output files, but do not submit anything. For testing'
            print '              -s or --skip-pythia       skip running pythia'
            print '              -h or --help              print this again'
            print '  Enjoy!'
            sys.exit(0)
   
      
      print 'doPythia = ' + str(self.doPythia)
      print 'submit   = ' + str(self.submit)




      # Assigns a unique task number.  This is a sixteen digit number.

      # For non-PanDA submission the top ten digits represent
      # the time of submission, with the bottom six being random
      # For PanDA submission, the top ten digits are the PanDA
      # job ID and the next six are zeroes.
  
 
      pandaID = os.environ.get('PandaID', None)
  
      if ( pandaID != None):
         taskID = pandaID + '000000'
      else:
         taskID=str(int(time.time()*1000000));

      #    Find out who is submitting this job and place it in 
      #   the variable 'user'. For some reason prodUserID
      #   is not set.

      self.user = os.environ.get('USER','nobody')
      if(self.user == 'apf'):
         self.user= os.environ.get('prodUserID','nobody')
      self.jobID = taskID + '0'

      # Converts the input strings that are numbers into numbers

      try :
         self.warmup_1 = int(warmup_1_string) 
      except ValueError:
         print 'Events for Warmup Phase 1 must be an integer'
         sys.exit(1)

      try :
         self.warmup_2 = int(warmup_2_string) 
      except ValueError:
         print 'Events for Warmup Phase 2 must be an integer'
         sys.exit(1)

      try :
         self.warmup_rounds = int(warmup_rounds_string) 
      except ValueError:
         print 'Rounds for Warmup must be an integer'
         sys.exit(1)

      try :
         self.n = int(event_number_string)
      except:
         print 'Number of events to generate is mandatory and must be an integer'
         sys.exit(1)

      # Checks to see that the process name is understood
  
      legal_process_names=['wqq', 'zqq', 'wcjet', 'wjet', 'zjet', 'hjet',
                           'vbjet', '2Q', '4Q', 'QQh', 'Njet', 'phjet', 
		                     'top', 'wphjet', 'wphqq', '2Qph' 
                          ]
      if not self.process_name in legal_process_names:
         print 'Unknown Alpgen process'
         sys.exit(1)
  
      return






class Shape:
# A place to store information on the various
#  job shapes.  It's convenient to do some
#  XML building here.


 def __init__(self, computer, nodes,  
             cores_used, events, sequence):


   if computer.name.find('vesta') >= 0:
      nodes = 32
      cores_used = 32
   elif computer.name.find('mira') >= 0:
      nodes = 512
      cores_used = 32
   else:
      print 'ERROR in Shape, nodes not an acceptable value\n'
      sys.exit(199)
   
   total_cores = nodes * cores_used
   n = int(math.ceil(events/total_cores))


   xml = '\t\t<shape>\n'
   xml = xml + '\t\t\t<nodes>' + str(nodes) + '</nodes>\n'
   xml = xml + '\t\t\t<cores_per_node>' + str(cores_used) 
   xml = xml + '</cores_per_node>\n'   

   if computer.architecture == 'bgp':
     if cores_used == 4:   
      xml = xml + '\t\t\t<flag>--mode=vn</flag>\n'
     if cores_used == 2:   
      xml = xml + '\t\t\t<flag>--mode=dual</flag>\n'
     if cores_used == 1:   
      xml = xml + '\t\t\t<flag>--mode=smp</flag>\n'

   elif computer.architecture == 'bgq':
     if   cores_used > 32:
      xml = xml + '\t\t\t<flag>--mode c64</flag>\n'
     elif cores_used > 16:
      xml = xml + '\t\t\t<flag>--mode c32</flag>\n'
     elif cores_used >  8:   
      xml = xml + '\t\t\t<flag>--mode c16</flag>\n'
     elif cores_used >  4:   
      xml = xml + '\t\t\t<flag>--mode c8</flag>\n'
     elif cores_used >  2:   
      xml = xml + '\t\t\t<flag>--mode c4</flag>\n'
     elif cores_used == 2:   
      xml = xml + '\t\t\t<flag>--mode c2</flag>\n'
     elif cores_used == 1:   
      xml = xml + '\t\t\t<flag>--mode c1</flag>\n'
     else:
      print 'Internal error'
      sys.exit(200)
   else:
     xml = xml + '\t\t\t<flag> </flag>\n'

   xml = xml + '\t\t\t<presubmit_script>presubmit.sh \n'
   xml = xml + '\t\t\t\t<presubmit_parameters>'
   xml = xml + str(total_cores) + ' ' 
   xml = xml + 'input.1.' + str(sequence)
   xml = xml + '</presubmit_parameters>\n'
   xml = xml + '\t\t\t</presubmit_script>\n'
   xml = xml + '\t\t\t<postsubmit_script>postsubmit.sh \n'
   xml = xml + '\t\t\t\t<postsubmit_parameters></postsubmit_parameters>\n'
   xml = xml + '\t\t\t</postsubmit_script>\n'
   xml = xml + '\t\t</shape>\n'
   
   self.xml = xml
   self.machine = computer
   self.seq = sequence
   self.n = n
   self.cores = total_cores

   if init.verbose:
     print xml


class ProcessInputFile:
# Parses the input file; the file named 'input'
#  and builds the input files for the three
#  separate phases

 def __init__(self, grid1, grid2, rounds):

   try:
     input = open('input','r')
   except:
     print 'Could not open input file, aborting'
     sys.exit(2)

   if os.path.exists(rundir+'input.0'):
     print 'File input.0 exists and will be overwritten.'
   try:
     output0 = open(rundir+'input.0','w')
   except:
     print 'Could not open input.0 file, aborting'
     sys.exit(2)

   if os.path.exists(rundir+'input.1.baseline'):
     print 'File input.1.baseline exists and will be overwritten.'
   try:
     output1= open(rundir+'input.1.baseline','w')
   except:
     print 'Could not open input.1.baseline file, aborting'
     sys.exit(2)

   if os.path.exists(rundir+'input.2'):
     print 'File input.2 exists and will be overwritten.'
   try:
     output2= open(rundir+'input.2','w')
   except:
     print 'Could not open input.2 file, aborting'
     sys.exit(2)     

# Line 1 - set mode
   try:
     input.readline()
   except:
     print 'Premature end of input file'
     sys.exit(2)  
   output0.write('0               ! imode\n')
   output1.write('1               ! imode\n')
   output2.write('2               ! imode\n')

# Line 2 - force name to alpout for internal use
   try:
      line=input.readline()
   except:
     print 'Premature end of input file'
     sys.exit(2)      
   line=line.strip()
   comment_position=line.find("!")
   if(comment_position > 0 ):
      self.label = line[:comment_position-1].rstrip()
   output0.write('alpout          ! label for files\n')   
   output1.write('alpout          ! label for files\n')   
   output2.write('alpout          ! label for files\n')         

# Line 3 - select the grid generation mode
#
# The fact that mode 2 has grid selection 1 is not a mistake

   try:
      line=input.readline()
   except:
     print 'Premature end of input file'
     sys.exit(2)  

   output0.write('0               !')   
   output0.write(' start with: 0=new grid, 1=previous warmup grid,')
   output0.write('  2=previous generation grid\n')

   output1.write('2               !') 
   output1.write(' start with: 0=new grid, 1=previous warmup grid,')
   output1.write('  2=previous generation grid\n')
   
   output2.write('1               !') 
   output2.write(' start with: 0=new grid, 1=previous warmup grid,')
   output2.write('  2=previous generation grid\n')
 
# Line 4 - set the round 1 warmup parameters 
 
   try:
      line=input.readline()
   except:
     print 'Premature end of input file'
     sys.exit(2)  
   output0.write('%d %d  ! Nevents/iteration,  N(warm-up iterations)\n'
                  %  (grid1, rounds))
   output1.write('0 0  ! Nevents/iteration,  N(warm-up iterations)\n')
   output2.write('0 0  ! Nevents/iteration,  N(warm-up iterations)\n')

# Line 5 - set the round 2 warmup parameters 
#   At this point, we don't know how many events to generate
#   per node.  The XXX's will need to be replaced later.

   try:
      line=input.readline()
   except:
     print 'Premature end of input file'
     sys.exit(2)  

   output0.write('%d   |   Nevents generated after warm-up\n' % grid2)
   output1.write('XXXXXXXXX   |   Nevents generated after warm-up\n')
   output2.write('999999999   |   Nevents generated after warm-up\n')
 
# Remaining lines
   self.ndns = -1
   while(len(line) > 0):
    line=input.readline()
    if(line[:5]=='njets'):
      comment_position=line.find('!')
      if(comment_position > 0):
        njets_text=line[6:comment_position-1].rstrip()
      else:
        njets_text=line[6:].rstrip()
      self.njets = int(njets_text)
      output0.write(line)
      output1.write(line)
      output2.write(line)
    elif(line[:4]=='ndns'):
      comment_position=line.find('!')
      if(comment_position > 0):
        ndns_text=line[5:comment_position-1].rstrip()
      else:
        ndns_text=line[5:].rstrip()
      self.ndns = int(ndns_text)
      output0.write(line)
      output1.write(line)
      output2.write(line)
    else:
      output0.write(line)
      output1.write(line)
      output2.write(line)

   input.close()  
   output0.close()
   output1.close()
   output2.close()


def runtime(events, arch, proc, njets):
# Estimates the run-time needed for a given job
#  Right now, this is largely fantasy, with 
#  lousy measurements for the BGP systems
#  and even worse extrapolations for the
#  other ones

  if arch=='bgp' :
   if proc == 'wjets':
    if njets == 0 :
      base_time = 300
      base_events = 400000    
    elif njets == 1:
      base_time = 507
      base_events = 150000
    elif njets == 3:
      base_time = 1236
      base_events = 300000
    elif njets == 4:
      base_time = 2725
      base_events = 2000000
    elif njets == 5:
      base_time = 3917
      base_events = 2500000
    else:
      base_time = 300
      base_events = 400000
   elif proc == 'zjets':
    if njets == 0 :
      base_time = 300
      base_events = 400000    
    elif njets == 1:
      base_time = 507
      base_events = 150000
    elif njets == 3:
      base_time = 1236
      base_events = 300000
    elif njets == 4:
      base_time = 2725
      base_events = 2000000
    elif njets == 5:
      base_time = 3917
      base_events = 2500000
    else:
      base_time = 300
      base_events = 400000
   else:
    if njets == 0 :
      base_time = 300
      base_events = 400000    
    elif njets == 1:
      base_time = 507
      base_events = 150000
    elif njets == 3:
      base_time = 1236
      base_events = 300000
    elif njets == 4:
      base_time = 2725
      base_events = 2000000
    elif njets == 5:
      base_time = 3917
      base_events = 2500000
    else:
      base_time = 300
      base_events = 400000
    seconds = events * base_time / base_events

  elif arch=='bgq' :
   if proc == 'wjets':
    if njets == 0 :
      base_time = 300
      base_events = 400000    
    elif njets == 1:
      base_time = 507
      base_events = 150000
    elif njets == 3:
      base_time = 1236
      base_events = 300000
    elif njets == 4:
      base_time = 2725
      base_events = 2000000
    elif njets == 5:
      base_time = 3917
      base_events = 2500000
    else:
      base_time = 300
      base_events = 400000
   elif proc == 'zjets':
    if njets == 0 :
      base_time = 300
      base_events = 400000    
    elif njets == 1:
      base_time = 507
      base_events = 150000
    elif njets == 3:
      base_time = 1236
      base_events = 300000
    elif njets == 4:
      base_time = 2725
      base_events = 2000000
    elif njets == 5:
      base_time = 3917
      base_events = 2500000
    else:
      base_time = 300
      base_events = 400000
   else:
    if njets == 0 :
      base_time = 300
      base_events = 400000    
    elif njets == 1:
      base_time = 507
      base_events = 150000
    elif njets == 3:
      base_time = 1236
      base_events = 300000
    elif njets == 4:
      base_time = 2725
      base_events = 2000000
    elif njets == 5:
      base_time = 3917
      base_events = 2500000
    else:
      base_time = 300
      base_events = 400000
    seconds = events * base_time / base_events / 2.0

  elif arch=='x86' :
   if proc == 'wjets':
    if njets == 0 :
      base_time = 300
      base_events = 400000    
    elif njets == 1:
      base_time = 507
      base_events = 150000
    elif njets == 3:
      base_time = 1236
      base_events = 300000
    elif njets == 4:
      base_time = 2725
      base_events = 2000000
    elif njets == 5:
      base_time = 3917
      base_events = 2500000
    else:
      base_time = 300
      base_events = 400000
   elif proc == 'zjets':
    if njets == 0 :
      base_time = 300
      base_events = 400000    
    elif njets == 1:
      base_time = 507
      base_events = 150000
    elif njets == 3:
      base_time = 1236
      base_events = 300000
    elif njets == 4:
      base_time = 2725
      base_events = 2000000
    elif njets == 5:
      base_time = 3917
      base_events = 2500000
    else:
      base_time = 300
      base_events = 400000
   else:
    if njets == 0 :
      base_time = 300
      base_events = 400000    
    elif njets == 1:
      base_time = 507
      base_events = 150000
    elif njets == 3:
      base_time = 1236
      base_events = 300000
    elif njets == 4:
      base_time = 2725
      base_events = 2000000
    elif njets == 5:
      base_time = 3917
      base_events = 2500000
    else:
      base_time = 300
      base_events = 400000
    seconds = events * base_time / base_events / 2.5
  else:
    print 'Internal error'
    sys.exit(200)  
      
  return(seconds)






##########################################################
#             Main program starts here                   #
########################################################## 

import stat, subprocess


init=Initialize();
if(init.verbose):
  print 'User name       :',init.user
  print 'Process name    :',init.process_name  
  print 'Job ID          :',init.jobID    
  print 'Warmup 1        :',init.warmup_1
  print 'Warmup 2        :',init.warmup_2
  print 'Warmup Rounds   :',init.warmup_rounds
  print 'Machines        :',init.machine_list
  print 'Events          :',init.n
  print 'Test Mode Flag  :',init.testmode


topdir=os.getcwd()+'/'
rundir=os.getcwd()+'/rundir/'+init.jobID+'/'

# Hard-coded values are here:

safety_factor = 1.5
rabbit_hutch = 'python /users/hpcusers/rabbitmq/'

# Hard-coded values for the HPC wait stage 
#   (mostly time-outs) are here:
# Defaults are: 
#   Poll every minute
#     for up to two hours
#   If the log file appears without a message
#   being received, after 5 minutes, conclude
#   the job has ended.

hpcwait_polltime= 60
hpcwait_limit = 86400
hpcwait_waittime = 7200

# Hard-coded directories

outdir='/grid/atlas/hpc/transfer/from_hpc/'+init.jobID+'/'
indir='/grid/atlas/hpc/transfer/to_hpc/'+init.jobID+'/'
wwwdir='/grid/www/hpc/jobs/'
pdfdir='/grid/atlas/hpc/common/alpgen/pdf/'
alpexedir='/grid/atlas/hpc/common/alpgen/executables/'
urlhead='http://atlasgridftp02.hep.anl.gov:40000/hpc/jobs/'
scriptdir='/grid/atlas/hpc/common/alpgen/scripts/'
pyexe='/grid/atlas/hpc/common/pythia/executables/runPythiaOnAlpgen'

# And temporarily overwritten ones go here

# Match them with what was requested
requested_computers = []
for s in init.machine_list:
   try:
      requested_computers.append( computers[s])
   except IndexError,e:
      print 'computer not found, ' + s + ' exception: ' + str(e)
      sys.exit(4444)


# Make and populate the directories

if not os.path.exists(rundir) :
  os.makedirs(rundir)   

if not os.path.exists(indir):
  os.makedirs(indir)   
  
if not os.path.exists(outdir):
  os.makedirs(outdir)   

if os.path.exists(alpexedir+init.process_name+'gen'):
  shutil.copy(alpexedir+init.process_name+'gen',rundir);
else:
  print "Can't find Alpgen executable"
  sys.exit(1)

if os.path.exists(pyexe):
   shutil.copy(pyexe,rundir)
else:
   print "Can't find Pythia executable"
   sys.exit(1)

# TEMPORARY REMOVE!!!
#shutil.copy('/users/hpcusers/generators/alpgen-v214/zjetwork/zjetgen',rundir)


# Get the input file
infile=ProcessInputFile(init.warmup_1, 
                        init.warmup_2, 
			init.warmup_rounds)

# Build the PDF dictionary
#  Since there's only one of them, might as well have
#  it in the main body of the code
#
# ndns=101 is broken in Alpgen, but is included anyway

pdfdict = {  1:'cteq4m.tbl', 
             2:'cteq4l.tbl',
             3:'cteq4hj.tbl',
             4:'cteq5m.tbl',	    	    
             5:'cteq5l.tbl',
             6:'cteq5hj.tbl',
             7:'cteq6m.tbl',
             8:'cteq6l.tbl',
             9:'cteq6l1.tbl',	    
	    10:'cteq61.00.tbl',
	    11:'cteq61.01.tbl',	    
	    12:'cteq61.02.tbl',
	    13:'cteq61.03.tbl',	  	    
	    14:'cteq61.04.tbl',
	    15:'cteq61.05.tbl',	    
	    16:'cteq61.06.tbl',
	    17:'cteq61.07.tbl',	
	    18:'cteq61.08.tbl',
	    19:'cteq61.09.tbl',		    
	    20:'cteq61.10.tbl',
	    21:'cteq61.11.tbl',	    
	    22:'cteq61.12.tbl',
	    23:'cteq61.13.tbl',	  	    
	    24:'cteq61.14.tbl',
	    25:'cteq61.15.tbl',	    
	    26:'cteq61.16.tbl',
	    27:'cteq61.17.tbl',	
	    28:'cteq61.18.tbl',
	    29:'cteq61.19.tbl',	    
	    30:'cteq61.20.tbl',
	    31:'cteq61.21.tbl',	    
	    32:'cteq61.22.tbl',
	    33:'cteq61.23.tbl',	  	    
	    34:'cteq61.24.tbl',
	    35:'cteq61.25.tbl',	    
	    36:'cteq61.26.tbl',
	    37:'cteq61.27.tbl',	
	    38:'cteq61.28.tbl',
	    39:'cteq61.29.tbl',	
	    40:'cteq61.30.tbl',
	    41:'cteq61.31.tbl',	    
	    42:'cteq61.32.tbl',
	    43:'cteq61.33.tbl',	  	    
	    44:'cteq61.34.tbl',
	    45:'cteq61.35.tbl',	    
	    46:'cteq61.36.tbl',
	    47:'cteq61.37.tbl',	
	    48:'cteq61.38.tbl',
	    49:'cteq61.39.tbl',
	   101:'cor01.dat',	
           102:'alf119.dat',	    
           103:'alf117.dat',
           104:'alf121.dat',	    
           105:'j121.dat',	   
           106:'lo2002.dat' }	   


# Get the PDF file
if(infile.ndns > 0):
 try:
   pdfname = pdfdict[infile.ndns]
 except:
   print 'Unknown PDF set'
   sys.exit(2)
 if os.path.exists(pdfdir+pdfname):
   shutil.copy(pdfdir+pdfname, rundir)
 else:
   print 'Unknown PDF set'
   sys.exit(2)



# Start building the shapes 
shapes = []
seq = 1
for machine in requested_computers:
   nodes = machine.min_nodes
   shapes.append(
                 Shape(machine, 
                       nodes, 
	                    machine.cores_per_node,
                       init.n,
		                 seq)
                )



# Build the input files, one per shape
#   Only the fifth line (number of events)
#    is different

try:
  input = open(rundir+'input.1.baseline','r')
except:
  print 'Could not open input.baseline file, aborting'
  sys.exit(3)
for s in shapes:
  try:
    output = open(rundir+'input.1.'+str(s.seq),'w')
  except:
    print 'Could not open input.1 file for writing, aborting'
    sys.exit(3)
  output.write(input.readline())
  output.write(input.readline())
  output.write(input.readline())
  output.write(input.readline()) 
  input.readline()
  output.write(str(s.n)
       +'       |   Nevents generated after warm-up\n')  
  line=input.readline()
  while(len(line) > 0):
    output.write(line)   
    line=input.readline()
  output.close()
  input.seek(0)
input.close()

# If this is testmode run, make an input.1

if init.testmode:
  for s in shapes:
     if s.machine.name == 'frontend':
        shutil.copy(rundir+"input.1."+str(s.seq),rundir+"input.1")

# Build the pre-submit scripts
#  Disabled 13-Jun-2013

#for s in shapes:
#  filename = rundir+'prep_'+str(s.seq)+'.sh'
#  try:
#    output = open(filename,'w')
#  except:
#    print 'Could not open prep file for writing, aborting'
#    sys.exit(3)
#  output.write('#! /bin/bash\n')
#  output.write('cp input.1.'+str(s.seq)+' input.1\n')
#  output.write('mkdir work\n')
#  output.write('for ((i=0; i<%d; i++ ))\n' % s.cores)
#  output.write('do\n')
#  output.write(" istr=`printf ")
#  output.write('"%05d" $i')
#  output.write(" `\n")
#  output.write(' mkdir work/$istr \n')
#  output.write(' cp alpout.grid1 work/$istr/alpout.grid1 \n')
#  output.write(' cp alpout.grid2 work/$istr/alpout.grid2 \n')
#  output.write(' cp *.tbl  work/$istr/ \n')
#  output.write('done\n')
#  output.close()
#  os.chmod(filename, 
#           stat.S_IRWXU | stat.S_IRWXG | 
#           stat.S_IROTH | stat.S_IXOTH)

# Build the XML file

filename = rundir+'alpgen.xml'
try:
  xml = open(filename,'w')
except:
  print 'Could not open XML file for writing, aborting'
  sys.exit(3)
xml.write('<?xml version=\"1.0\" encoding=\"UTF-8\" ?> \n')
xml.write('<HPCjob>\n')
xml.write('\t<program>Alpgen\n')
xml.write('\t\t<process>'+init.process_name+'</process>\n')
xml.write('\t</program>\n\n')
xml.write('\t<identity>\n')
xml.write('\t\t'+init.user+'\n')
xml.write('\t</identity>\n')
for r in requested_computers:
  xml.write('\t<computer>\n')
  xml.write('\t\t<name>'+r.name+'</name>\n')
  jobtime = int(math.ceil(safety_factor
            *runtime(init.n, r.architecture, 
	             init.process_name, infile.njets)/3600))
  if r.name.find('vesta') >= 0:
   jobtime = 1 # hours
  elif r.name.find('mira') >= 0:
   jobtime = 1 # hours
  else:
   print 'ERROR finding machine name.\n'
   sys.exit(201)
  xml.write('\t\t<cpu_hours>'+str(jobtime)+'</cpu_hours>\n')
  for s in shapes:
    if s.machine == r:
      xml.write(s.xml)   
  xml.write('\t</computer>\n')   
xml.write('\t<transfer>\n')
xml.write('\t\t<inputfile>\n')
xml.write('\t\t\tgsiftp://atlasgridftp02.hep.anl.gov/'+indir+'\n')    		
xml.write('\t\t</inputfile>\n')
xml.write('\t\t<outputdir>\n')
xml.write('\t\t\tgsiftp://atlasgridftp02.hep.anl.gov/'+outdir+'\n') 
xml.write('\t\t</outputdir>\n')
xml.write('\t</transfer>\n')
xml.write('</HPCjob>\n')
xml.close()
	 

# Now write the condor scripts and DAG
print 'write condor cmd files for phase0/1/2/3'
phase0 = CondorJob(
                    comment               = 'run Alpgen warmup',
                    executable            = init.process_name+'gen',
                    input                 = 'input.0',
                    output                = 'warmup.out',
                    log                   = 'warmup.log',
                    error                 = 'warmup.err',
                    transfer_input_files  = pdfname,
                    transfer_output_files = 'alpout.grid1,alpout.grid2'
                   )
phase0.write(rundir+'warmup.cmd')

phase1 = CondorJob(
                    comment               = 'run hpcwait',
                    executable            = 'hpcwait.sh',
                    arguments             = str(init.jobID) + '  ' + outdir,
                    output                = 'hpcwait.out',
                    log                   = 'hpcwait.log',
                    error                 = 'hpcwait.err',
                   )
phase1.write(rundir+'hpcwait.cmd')

phase2 = CondorJob(
                    comment               = 'run Alpgen unweighting',
                    executable            = init.process_name+'gen',
                    input                 = 'input.2',
                    output                = 'unweight.out',
                    log                   = 'unweight.log',
                    error                 = 'unweight.err',
                    transfer_input_files  = pdfname +',alpout.grid1, alpout.grid2, alpout.wgt, alpout.par',
                    transfer_output_files = 'alpout.unw,alpout_unw.par',
                   )
phase2.write(rundir+'unweight.cmd')

phase3 = CondorJob(
                    comment               = 'run Pythia hadronization on the Alpgen Unweighted events',
                    executable            = 'runPythiaOnAlpgen',
                    output                = 'pythia.out',
                    log                   = 'pythia.log',
                    error                 = 'pythia.err',
                    transfer_input_files  = 'alpout.unw,alpout_unw.par,pythia.cmnd',
                    transfer_output_files = 'pythia_output.hepmc',
                   )
phase3.write(rundir+'pythia.cmd')

print 'writing alpgen.dag file'
try:
   dag = open(rundir+'alpgen.dag','w')
except IOError,e:
   print 'ERROR opening file alpgen.dag for writing, exception: ' + str(e)
   sys.exit(3)

dag.write('Job A warmup.cmd\n')
dag.write('Script Post A place_input_files.sh\n')
dag.write('Job B hpcwait.cmd\n')
dag.write('Job C combine.cmd\n')
dag.write('Job D unweight.cmd\n')
if init.doPythia:
   dag.write('Job E pythia.cmd\n')

dag.write('Parent A Child B')
dag.write(' C D')
if init.doPythia:
   dag.write(' E')
dag.write('\n')
dag.write('Parent B Child C D')
if init.doPythia:
   dag.write(' E')
dag.write('\n')
dag.write('Parent C Child D')
if init.doPythia:
   dag.write(' E')
dag.write('\n')
if init.doPythia:
   dag.write('Parent D Child E\n')

dag.close()

# Get additional executables
print 'copy scripts to run directory: ',rundir
shutil.copy(scriptdir+'alpcounter', rundir)
shutil.copy(scriptdir+'presubmit.sh', rundir)
shutil.copy(scriptdir+'postsubmit.sh',rundir)

print 'writing hpcwait script'
try:
  hpcwait = open(rundir+'hpcwait.sh','w')
except:
  print 'Could not open HPC waiting script for writing, aborting'
  sys.exit(3)

hpcwait.write('#! /bin/bash\n\n')
hpcwait.write(
     '# Waits for HPC job to finish by waiting for the log file to appear\n')
hpcwait.write(
      '#  This assumes that the log file is the last to make it\n')
hpcwait.write(
      '# Arguments: the Job ID and the directory that should be polled \n\n')

hpcwait.write('POLLTIME='+str(hpcwait_polltime)+'\n')
hpcwait.write('LIMIT='+str(hpcwait_limit)+'\n')
hpcwait.write('WAITTIME='+str(hpcwait_waittime)+'\n\n')
hpcwait.write('JOBDONE=0\n')
hpcwait.write('LOGEXIST=0\n')
hpcwait.write('TIMER=1\n\n')

hpcwait.write('while [[ $TIMER -lt $LIMIT  &&  $JOBDONE -eq 0  ]]\n')
hpcwait.write('do\n')
for r in requested_computers:
  hpcwait.write(rabbit_hutch+
                'receive_message_on_hpc.py '+
                str(init.jobID)+'.finishedjob\n')
  hpcwait.write('  if [ $? == 0 ]\n')
  hpcwait.write('  then\n')
  hpcwait.write('     echo "Job Done"\n')
  hpcwait.write('     JOBDONE=1\n')
  hpcwait.write('  fi\n')
hpcwait.write('\n')

# What to do when we have indications of a completed job, but no message

hpcwait.write('  if [ -e $2/*.cobaltlog ]\n')
hpcwait.write('  then\n')
hpcwait.write('      echo "Log exists"\n')
hpcwait.write('      LOGEXIST=1\n')
hpcwait.write('     if [ $WAITTIME -le 0 ]\n')
hpcwait.write('     then\n')
hpcwait.write('          echo "Messages never arrived"\n')
hpcwait.write('          JOBDONE=1\n')
hpcwait.write('     fi\n')
hpcwait.write('     let WAITTIME=WAITTIME-POLLTIME\n')
hpcwait.write('  fi\n\n')

hpcwait.write('  if [ $TIMER -gt $LIMIT ]\n')
hpcwait.write('     then\n')
hpcwait.write('     echo "Timeout waiting for HPC"\n')
hpcwait.write('     '+rabbit_hutch+'delete_taskID_queue.py '+str(init.jobID)+'\n')
hpcwait.write('     exit 1;\n')
hpcwait.write('  fi\n')
hpcwait.write('  if [ $JOBDONE -eq 0 ]\n')
hpcwait.write('  then\n')
hpcwait.write('     sleep $POLLTIME\n')
hpcwait.write('     let TIMER=TIMER+POLLTIME\n')
hpcwait.write('  fi\n')
hpcwait.write('  echo $TIMER\n')

hpcwait.write('done\n\n')

hpcwait.write('  '+rabbit_hutch+'delete_taskID_queue.py '+str(init.jobID)+'\n')
hpcwait.write('  exit 0;\n')
hpcwait.close()
os.chmod(rundir+'hpcwait.sh', 
           stat.S_IRWXU | stat.S_IRWXG | 
           stat.S_IROTH | stat.S_IXOTH)






# Build placement script
print ' write place input files script'
try:
  pif = open(rundir+'place_input_files.sh','w')
except:
  print 'Could not open placement file for writing, aborting'
  sys.exit(3)
pif.write('#! /bin/bash \n')
pif.write('cp alpout.grid1 '+indir+'\n')
pif.write('cp alpout.grid2 '+indir+'\n')
pif.write('cp '+pdfdir+pdfname+' '+indir+'\n')
pif.write('cp presubmit.sh '+indir+'\n')
pif.write('chmod a+x '+indir+'/presubmit.sh\n')
pif.write('cp postsubmit.sh '+indir+'\n')
pif.write('chmod a+x '+indir+'/postsubmit.sh\n')
pif.write('cp input.* '+indir+'\n')
pif.write(rabbit_hutch+'create_taskID_queue.py '+init.jobID+'\n')

#  Since we don't have a control tower
#  we need some way to select a computer to run this on.  We 
#  use the last one in the list.

# Getjobs replaced with submit

if init.submit:
   pif.write(rabbit_hutch+'send_message.py '
            +requested_computers[-1].name+'.'
            +init.jobID+'.submit '+ urlhead + init.jobID + '.xml \n')

pif.close()
os.chmod(rundir+'place_input_files.sh', 
              stat.S_IRWXU | stat.S_IRWXG | 
              stat.S_IROTH | stat.S_IXOTH)


# Build combiner script and submit file
print 'write combine script'
try:
  combine = open(rundir+'combine.sh','w')
except:
  print 'Could not open combine script file for writing, aborting'
  sys.exit(3)
combine.write('#! /bin/bash\n')
combine.write('cp '+outdir+'alpout.wgt alpout.wgt\n')
combine.write('if [ -f '+outdir+'alpout.par ]\n')
combine.write('    then\n')
combine.write('    pwd \n')
combine.write('    ls -l \n')
combine.write('    cp '+outdir+'alpout.par alpout.par\n')
combine.write('fi\n')
combine.write('if [ -f '+outdir+'work/00000/alpout.par ]\n')
combine.write('    then\n')
combine.write('    cp '+outdir+'work/00000/alpout.par alpout.par\n')
combine.write('fi\n')
combine.write('./alpcounter \n')
combine.write('rm alpout.par\n')
combine.write('cp alpout.par.new alpout.par\n')
combine.close()
os.chmod(rundir+'combine.sh', 
           stat.S_IRWXU | stat.S_IRWXG | 
           stat.S_IROTH | stat.S_IXOTH)

print 'write combine condor file'
condor_combine = CondorJob(
                       comment               = 'run combine script',
                       executable            = 'combine.sh',
                       output                = 'combine.out',
                       log                   = 'combine.log',
                       error                 = 'combine.err',
                       transfer_input_files  = 'alpcounter',
                       transfer_output_files = 'alpout.wgt,alpout.par',
                      )
condor_combine.write(rundir+'combine.cmd')

# Build HPC fake submit file

if init.testmode:

  try:
    fake = open(rundir+'hpcfake.cmd','w')
  except:
    print 'Could not open HPC fake HTCondor file for writing, aborting'
    sys.exit(3)

  fake.write('####################\n')
  fake.write('universe      = vanilla\n')
  fake.write('executable    = '+init.process_name+'gen\n')
  fake.write('input         = input.1\n')
  fake.write('output        = hpcfake.out\n')
  fake.write('error         = hpcfake.err\n')
  fake.write('log           = hpcfake.log\n\n')
  fake.write('should_transfer_files = YES\n')
  fake.write('when_to_transfer_output = ON_EXIT\n')
  fake.write('transfer_input_files =  alpout.grid1, alpout.grid2')
  if(infile.ndns > 0):
     fake.write(','+pdfname)
  fake.write('\n') 
  fake.write('transfer_output_files = alpout.wgt, alpout.par\n')
  fake.write('transfer_output_remaps = \" alpout.wgt = ')
  fake.write(outdir+'alpout.00000.wgt ')
  fake.write('; alpout.par = ')
  fake.write(outdir+'alpout.par "\n')
  fake.write('queue\n')
  fake.write('####################\n')
  fake.close()

  try:
    fakedone = open(rundir+'hpcfakedone.sh','w')
  except:
    print 'Could not open HPC fake done script for writing, aborting'
    sys.exit(3)
  fakedone.write('#! /bin/bash\n')
  fakedone.write(rabbit_hutch+'send_message.py frontend.'+
            init.jobID+'.finishedjob '+
            urlhead+init.jobID+'.xml \n')
  fakedone.write('touch '+outdir+'alpgen.cobaltlog\n')
  fakedone.close()
  os.chmod(rundir+'hpcfakedone.sh', 
           stat.S_IRWXU | stat.S_IRWXG | 
           stat.S_IROTH | stat.S_IXOTH)  


# create pythia config file
if init.doPythia:
   print 'write pythia config file'
   try:
      pythiaconf = open(rundir+'pythia.cmnd','w')
   except:
     print 'Could not open pythia config file for writing, aborting'
     sys.exit(3)

   pythiaconf.write('! 1) Settings used in the main program.\n')
   pythiaconf.write('Main:numberOfEvents   = -1         ! number of events to generate (-1 for all)\n')
   pythiaconf.write('Main:timesAllowErrors = 3          ! how many aborts before run stops\n')
   pythiaconf.write('Main:spareMode1 = 0                ! skip n events at beginning of file\n')
   pythiaconf.write('\n')
   pythiaconf.write('! 2) Settings related to output in init(), next() and stat().\n')
   pythiaconf.write('Init:showChangedSettings = on      ! list changed settings\n') 
   pythiaconf.write('Init:showChangedParticleData = on  ! list changed particle data\n')
   pythiaconf.write('Next:numberCount       = 10        ! print message every n events\n')
   pythiaconf.write('Next:numberShowInfo    = 1         ! print event information n times\n')
   pythiaconf.write('Next:numberShowProcess = 1         ! print process record n times\n')
   pythiaconf.write('Next:numberShowEvent   = 1         ! print event record n times\n')
   pythiaconf.write('\n')
   pythiaconf.write('! 3) Enable matching\n')
   pythiaconf.write('JetMatching:merge = on\n')
   pythiaconf.write('\n')
   pythiaconf.write('! Alpgen run\n')
   pythiaconf.write('Alpgen:file = alpout\n')
   pythiaconf.write('Alpgen:setMLM = on\n')
   pythiaconf.write('JetMatching:scheme = 2\n')
   pythiaconf.write('JetMatching:exclusive = 1\n')
   pythiaconf.write('\n')
   pythiaconf.close()
   os.chmod(rundir+'pythia.cmnd',
            stat.S_IRWXU | stat.S_IRWXG |
            stat.S_IROTH | stat.S_IXOTH)

# create pythia step condor 


# Publish job parameters 
print 'copy xml file to webserver'
shutil.copy(rundir+'alpgen.xml', wwwdir+init.jobID+'.xml')

# Submit the jobs - we're willing to wait for the HPC to time out,
#   plus one hour
if init.submit:
   print 'submit dag job to condor'
   os.chdir( rundir )
   subprocess.call(['condor_submit_dag','alpgen.dag'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   print 'call condor_wait'
   subprocess.call(['condor_wait',
      '-wait',str(hpcwait_limit+3600),'alpgen.dag.dagman.log'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)


# This would be a good place to put some error checking


# Cleanup

#os.system('cp '+outdir+'*.cobaltlog ' + rundir)
if init.submit:
   print 'tar output'
   os.system('tar cvzf '+topdir+'logs.tgz'+
             ' *.log *.err *.out ')
   os.chdir(topdir)
   shutil.move(rundir+'alpout.unw',infile.label+'.unw')
   if init.doPythia:
      shutil.move(rundir+'pythia_output.hepmc','pythia_output.hepmc')
   os.remove(wwwdir+init.jobID+'.xml')
   #shutil.rmtree(rundir)
   #shutil.rmtree(indir)
   #shutil.rmtree(outdir)

#Do we want to do this?  I think we don't.
#os.remove('input')

sys.exit(0)



