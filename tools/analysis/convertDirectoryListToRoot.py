#!/usr/bin/env python
import ROOT
ROOT.gROOT.SetBatch(True)

import optparse,datetime,sys,time,array

# I haven't been able to get the filename to properly appear in the TTree. not sure what is going on there.

def main():
   parser = optparse.OptionParser(description='convert directory dump to TTree in a TFile')
   parser.add_option('-i','--input-file',dest='inputFilename',help='input directory list text file')
   parser.add_option('-o','--output-file',dest='outputFilename',help='output root filename',default='output.root')
   parser.add_option('-c','--cobalt-log',dest='cobaltlog',help='can provide the cobalt log file and use the job start time')
   parser.add_option('-a','--offset-hours',dest='offset_hours',help='if you need to add an offset to the hours',type='int')
   parser.add_option('-n','--new-ls',dest='new_ls',action='store_true',help='use the new verbose ls',default=False)
   options,args = parser.parse_args()

   if options.inputFilename is None:
      parser.error('Must specify an input file using -i')

   start_time = datetime.datetime(2014,1,1)
   offset = datetime.timedelta(hours=0)
   if options.offset_hours is not None:
      offset = datetime.timedelta(hours=options.offset_hours)
   # if cobalt log file given, set start time
   if options.cobaltlog is not None:
      cobalt_file = open(options.cobaltlog)
      for line in cobalt_file:
         if line.find('Initiating boot at') >= 0:
            start_time = datetime.datetime.strptime(line[0:24], '%a %b %d %H:%M:%S %Y')
            start_time = start_time - datetime.timedelta(hours=5) + offset
            break
   print 'start time = ' + str(start_time)
   
   infile = open(options.inputFilename)
   
   bytes          = array.array('L',[0])
   dateInSeconds  = array.array('L',[0])
   filename       = array.array('c','filename\0')

   # create TTree
   tree = ROOT.TTree('directoryList','directoryList')
   tree.Branch('bytes',bytes,'bytes/i')
   tree.Branch('dateInSeconds',dateInSeconds,'dateInSeconds/i')
   tree.Branch('filename',filename,'filename/C')
   
   secondsPerDay = 60*60*24

   for line in infile:
      parts = line.split()
      # if there are not 9 parts to the list then this is not a file or directory listing
      if len(parts) != 9:
         continue
      
      # identify the parts
      permissions = parts[0] # rwx------
      unknownA    = parts[1] # 1
      owner       = parts[2]
      group       = parts[3]
      bytes[0]    = int(parts[4])
      file        = parts[8]
      if options.new_ls:
         date        = parts[5]
         time        = parts[6]
         timezone    = parts[7]
         #dateTime = datetime.datetime.strptime(date + ' ' + time + ' ' + timezone,"%Y-%m-%d %H:%M:%S.%f %z")
         dateTime = datetime.datetime.strptime(date + ' ' + time.split('.')[0],"%Y-%m-%d %H:%M:%S")
      else:
         date = '%04d %s %02d %s' % (2014,parts[5],int(parts[6]),parts[7])
         dateTime = datetime.datetime.strptime(date, "%Y %b %d %H:%M")

      if ('.wgt' in file or '.unw' in file) and 'alpout.unw' not in file:
         #print str(dt)
         tdelta = (dateTime - start_time)
         dateInSeconds[0] = int(tdelta.days*secondsPerDay + tdelta.seconds)
         #print str(dateInSeconds[0]) + ' ' + str(dateTime)

         filename = array.array('c',file+'\0')
         tree.SetBranchAddress('filename',filename)

         tree.Fill()
      
         #print ("%010d %010d %s" % (bytes[0],dateInSeconds[0],filename)) 

   file = ROOT.TFile(options.outputFilename,'RECREATE')
   tree.Write(tree.GetName(),ROOT.TObject.kOverwrite)
   file.Close()

   return 0



if __name__ == '__main__':
   main()

