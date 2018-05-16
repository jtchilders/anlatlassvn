#!/usr/bin/env python
import os,sys,optparse,logging
logger = logging.getLogger(__name__)

import ROOT
ROOT.gROOT.SetBatch(True)
from AtlasStyle import SetAtlasStyle
SetAtlasStyle()

from root.CanvasFile import CanvasFile
from root.TGraph import TGraph
from sherpa.SherpaLogFile import SherpaLogFile

def main():
   logging.basicConfig(level=logging.INFO)

   parser = optparse.OptionParser(description='Use the log files from sherpa runs to plot run time data.')
   parser.add_option('-i','--input',dest='input',help='comma separated list of Sherpa log files.')
   parser.add_option('-o','--output',dest='output',help='Output filename. Extension determines file output: .ps .png .jpg',default=sys.argv[0] + '.ps')
   options,args = parser.parse_args()

   if options.input is None:
      parser.error('Must specify -i')

   log_filenames = options.input.split(',')
   if len(log_filenames) < 1:
      parser.error('syntax of input list is incorrect')

   logger.info(' Looping over these files: ' + str(log_filenames))
   logger.info(' ouput results to file: ' + options.output)

   g_init_time_vs_rank = TGraph('g_init_time_vs_rank')
   h_init_time_vs_rank = ROOT.TH1D('h_init_time_vs_rank',';number of ranks;time (s)',20,0,20)
   h_init_time_vs_rank.SetFillColor(ROOT.kRed)
   h_init_time_vs_rank.SetLineWidth(0)
   g_init_time_vs_thread = TGraph('g_init_time_vs_thread')
   h_init_time_vs_thread = ROOT.TH1D('h_init_time_vs_thread',';number of threads;time (s)',20,0,20)
   h_init_time_vs_thread.SetFillColor(ROOT.kRed)
   h_init_time_vs_thread.SetLineWidth(0)

   g_evtgen_time_vs_rank = TGraph('g_evtgen_time_vs_rank')
   h_evtgen_time_vs_rank = ROOT.TH1D('h_evtgen_time_vs_rank',';number of ranks;time (s)',20,0,20)
   h_evtgen_time_vs_rank.SetFillColor(ROOT.kBlue)
   h_evtgen_time_vs_rank.SetLineWidth(0)
   g_evtgen_time_vs_thread = TGraph('g_evtgen_time_vs_thread')
   h_evtgen_time_vs_thread = ROOT.TH1D('h_evtgen_time_vs_thread',';number of threads;time (s)',20,0,20)
   h_evtgen_time_vs_thread.SetFillColor(ROOT.kBlue)
   h_evtgen_time_vs_thread.SetLineWidth(0)

   h_evtgen_evt_rate = ROOT.TH1D('h_evtgen_evt_rate',';number of ranks; events per second per rank',20,0,20)

   for log_filename in log_filenames:
      try:
         log = SherpaLogFile.read_file(log_filename)
      except:
         logger.exception('Error reading log file ' + log_filename)
         raise

      g_init_time_vs_rank.AddPoint(log.num_ranks,log.init_sec)
      h_init_time_vs_rank.Fill(log.num_ranks,log.init_sec)
      g_init_time_vs_thread.AddPoint(log.num_threads,log.init_sec)
      h_init_time_vs_thread.Fill(log.num_threads,log.init_sec)
      g_evtgen_time_vs_rank.AddPoint(log.num_ranks,log.evt_gen_sec)
      h_evtgen_time_vs_rank.Fill(log.num_ranks,log.evt_gen_sec)
      g_evtgen_time_vs_thread.AddPoint(log.num_threads,log.evt_gen_sec)
      h_evtgen_time_vs_thread.Fill(log.num_threads,log.evt_gen_sec)

      h_evtgen_evt_rate.Fill(log.num_ranks*log.num_threads,log.num_events_per_rank*1./log.evt_gen_sec)

   
   cf = CanvasFile()
   cf.AddPage()
   cf.AddPlot(g_init_time_vs_rank.graph,'ap')
   cf.AddPage()
   cf.AddPlot(g_init_time_vs_thread.graph,'ap')
   cf.AddPage()
   cf.AddPlot(g_evtgen_time_vs_rank.graph,'ap')
   cf.AddPage()
   cf.AddPlot(g_evtgen_time_vs_thread.graph,'ap')

   g_init_time_vs_rank.SetAxisLabels('number of ranks','initialization time (s)')
   g_init_time_vs_thread.SetAxisLabels('number of threads','initialization time (s)')
   g_evtgen_time_vs_rank.SetAxisLabels('number of ranks','event generation time (s)')
   g_evtgen_time_vs_thread.SetAxisLabels('number of threads','event generation time (s)')

   cf.AddPage()
   stack_ranks = ROOT.THStack('stack_ranks',';number of ranks;time (s)')
   stack_ranks.Add(h_init_time_vs_rank)
   stack_ranks.Add(h_evtgen_time_vs_rank)
   cf.AddPlot(stack_ranks)

   cf.AddPage()
   stack_threads = ROOT.THStack('stack_threads',';number of threads;time (s)')
   stack_threads.Add(h_init_time_vs_thread)
   stack_threads.Add(h_evtgen_time_vs_thread)
   cf.AddPlot(stack_threads)

   cf.AddPage()
   cf.AddPlot(h_evtgen_evt_rate)
   

   cf.SaveAs(options.output)


   


if __name__ == "__main__":
   main()
