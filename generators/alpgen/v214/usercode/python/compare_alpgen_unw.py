#!/usr/bin/env python
import ROOT
ROOT.gROOT.SetBatch(True)
import sys,logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
from CanvasFile import CanvasFile
from AtlasStyle import SetAtlasStyle
from RatioCanvas import RatioCanvas

SetAtlasStyle()


filenames = [
             '/grid/atlas/hpc/argo/jobs/14042388640693090/alpout.root',
             '/grid/atlas/hpc/argo/jobs/14042422704408350/alpout.root',
             '/grid/atlas/hpc/argo/jobs/14036156039309510/alpout.root',
            ]
legend_labels = [
                 '4096x32',
                 '32x32',
                 '4096x32',
                ]

normalize = True

ratio_plots_filename = 'ratios.ps'


files = []
for filename in filenames:
   try:
      files.append(ROOT.TFile.Open(filename))
   except e:
      print 'ERROR opening',filename,', exception:',str(e)
      sys.exit(-1)

canvas = ROOT.TCanvas('canvas','canvas',0,0,800,600)
canvas.SetMargin(0,0,0.6,0.3)

for i in range(files[0].GetListOfKeys().GetSize()):
   plot_name = files[0].GetListOfKeys().At(i).GetName()
   print 'plotting',plot_name

   # set legend location
   legend_x1 = 0.7
   legend_x2 = 0.98
   legend_y1 = 0.9
   legend_y2 = 0.7
   min_y = 0
   log_y = False
   min_x = None
   max_x = None
   if plot_name.find('phi') >= 0:
      legend_y1 = 0.5
      legend_y2 = 0.3
   if plot_name.find('z_pt') >= 0:
      max_x = 80
      min_x = 25

   ratioCanvas = RatioCanvas(canvas=canvas,
                             ratioPad_y_ndiv = 204,
                             ratioPad_y_title = 'X/512',
                             ratioPad_bottomMargin = 0.5,
                             legend_chi2 = True,
                            )
   # get all the same named plots from each file
   for j in range(len(files)):
      plot = files[j].Get(plot_name)
      plot.SetDirectory(0)
      plot.Sumw2()
      plot.SetLineColor(ROOT.kBlack + j)
      plot.SetMarkerColor(ROOT.kBlack + j)
      plot.SetMarkerStyle(21+j)
      plot.SetMarkerSize(0.4)
      plot.SetLineWidth(len(files)-j)
      if normalize:
         plot.Scale(1.0/plot.Integral('width'))
      if j == 0:
         ratioCanvas.AddPlot(plot,'e',legend_labels[j],'lp',True)
      else:
         ratioCanvas.AddPlot(plot,'e',legend_labels[j],'lp')
   
   if i == 0:
      ratioCanvas.Draw(ratio_plots_filename+'(')
   elif i == files[0].GetListOfKeys().GetSize() - 1:
      ratioCanvas.Draw(ratio_plots_filename+')')
   else:
      ratioCanvas.Draw(ratio_plots_filename)

   
      
   

