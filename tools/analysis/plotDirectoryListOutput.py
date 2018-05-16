#!/usr/bin/env python
from ROOT import gROOT,TFile,TH1D,TH2D
gROOT.SetBatch(True)
import CanvasFile
import AtlasStyle
AtlasStyle.SetAtlasStyle()
import sys,os

def main():
   # pass a list of files on the command line

   if len(sys.argv) == 1:
      print 'usage: ' + sys.argv[0] + ' file1.root:legend_label1 [file2.root:legend_label2 ... ]'
      return
   
   textlist = sys.argv[1:]

   canvas = CanvasFile.CanvasFile()
   canvas.AddPage(legend_x1 = 0.2,legend_y1 = 0.9,legend_x2 = 0.5,legend_y2 = 0.6,auto_max_y = True,log_y = True)
   
   i = 0
   for text in textlist:
      split_text = text.split(':')
      filename = split_text[0]
      legend_label = None
      if len(split_text) > 1:
         legend_label = split_text[1]
      print ' plotting file ' + filename + ' with label ' + str(legend_label)
      file = TFile(filename)
      tree = file.Get('directoryList')
      entries = tree.GetEntriesFast()

      rank_file_times = TH1D('rank_file_times_'+str(i),';minutes from start of job',160,20,180)
      rank_file_times.SetDirectory(0)
      rank_file_times.SetLineColor(1+i)

      wgt_size_vs_time = TH2D('size_vs_time_'+str(i),
                              ';weighted event file size;minutes from start of job',
                              100,0,1500000,100,20,120)
      wgt_size_vs_time.SetDirectory(0)

      i += 1

      entry = 0
      while entry < entries:
         tree.GetEntry(entry)
         dateInMinutes = tree.dateInSeconds/60
         rank_file_times.Fill(dateInMinutes)
         if tree.filename.find('.wgt') >= 0:
            wgt_size_vs_time.Fill(tree.bytes,dateInMinutes)

         entry += 1

      canvas.AddPlot(rank_file_times,legend_label=legend_label,legend_options='l',page_number=0)
      canvas.AddPage()
      canvas.AddPlot(wgt_size_vs_time,'colz')

   canvas.SaveAs('rank_file_times.ps')



if __name__ == '__main__':
   main()

