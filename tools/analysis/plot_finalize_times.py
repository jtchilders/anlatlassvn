#!/usr/bin/env python

import os,sys,subprocess,optparse,datetime,glob
from ROOT import TGraph,TCanvas,TH1D,gROOT
gROOT.SetBatch(True)
import AtlasStyle
AtlasStyle.SetAtlasStyle()

START_GEN_EVT  = 'Generating weighted events' # Thu Oct 2 16:05:59 UTC 2014
INITIALIZE_MPI = 'INITIALIZE RANK' # 00005915 AT  46.70999908
FINALIZE_MPI   = 'FINALIZE RANK'   # 00005915 AT  46.70999908
START_UNW      = 'Unweighting events' # Thu Oct 2 16:05:59 UTC 2014
END_SCRIPT     = 'Done unweighting events' # Thu Oct 2 16:05:59 UTC 2014 

def main():
   parser = optparse.OptionParser(description='extract rank finalized times')
   parser.add_option('-i','--input-log',dest='log_filename',help='name of the cobalt output log file.')
   options,args = parser.parse_args()

   if options.log_filename is None:
      parser.error('Must specify an input log file using -i')

   
   print ' getting information from file: ' + str(options.log_filename)

   wgt_gen_start  = get_line_list(options.log_filename,START_GEN_EVT)
   unw_start      = get_line_list(options.log_filename,START_UNW)
   end_unw        = get_line_list(options.log_filename,END_SCRIPT)

   
   output_file = open('rank_time.txt','w')
   missed_file = open('rank_missed.txt','w')
   
   i=0
   diff = 0
   graph = TGraph()
   graph.GetXaxis().SetTitle('rank number')
   graph.GetYaxis().SetTitle('weighted evt gen finalize time')
   graph.SetMarkerStyle(1)
   histo = TH1D('histo',';weighted event generation rank finalize time',200,0,7000)
   for rank,seconds in rank_time.iteritems():
      if rank - i != diff:
         print ' rank: %10d' % rank
         missed_file.write('%10d' % rank)
      graph.SetPoint(i,rank,seconds)
      output_file.write("%10d%10d\n" % (rank,seconds))
      histo.Fill(seconds)
      diff = rank - i
      i+=1

   can = TCanvas("can","can",0,0,800,600)
   graph.Draw('ap')
   can.SaveAs('rank_time.ps(')
   can.SetLogy()
   histo.Draw()
   can.SaveAs('rank_time.ps)')

   print ' ranks read: ' + str(i)


def get_line_list(filename,string):
   
   try:
      line_list = []
      for line in open(filename):
         index = line.find(string)
         if index >= 0:
            line_list.append(line[0:-1])
            
      return line_list
            
            
   except:
      print ' exception received while trying to get rank time list: ' + str(sys.exc_info()[1])
      return None
   
if __name__ == '__main__':
   main()


