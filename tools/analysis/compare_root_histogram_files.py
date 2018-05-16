#!/usr/bin/env python
import os,sys,optparse,logging
logger = logging.getLogger(__name__)
# load up ROOT
try:
   import ROOT
except:
   raise Exception('ROOT not in PYTHONPATH')

from root.CanvasFile import CanvasFile

class RootFile:
   def __init__(self,filename):
      self.filename = filename
      file = ROOT.TFile(filename)
      self.th1 = {}
      for key in file.GetListOfKeys():
         obj = key.ReadObj()
         if isinstance(obj,ROOT.TH1):
            obj.SetDirectory(0)
            self.th1[obj.GetName()] = obj
      file.Close()
      file = None

class PageSettings:
   def __init__(self,page_number):
      self.page_number = page_number
      self.log_x = False
      self.log_y = False
      self.min_x = None
      self.max_x = None
      self.min_y = None
      self.max_y = None

def get_chi_square(histo):
   chisq = 0
   ndf = 0
   for bin_number in range(1,histo.GetNbinsX()+1):
      if(histo.GetBinError(bin_number) != 0.):
         chisq += (histo.GetBinContent(bin_number))**2 / (histo.GetBinError(bin_number))**2
         ndf += 1
   return chisq/ndf

def main():
   parser = optparse.OptionParser(description='Given a list of root files containing the same histograms, this will plot them together and test for agreement.')
   parser.add_option('-i','--input',dest='input',help='Comma separated list of input root files, each containing the same files.')
   parser.add_option('-o','--output',dest='output',help='Output TCanvas::SaveAs file name. The file postfix determines the format, for instance ".ps", ".jpg", or ".png" ',default='output.ps')
   parser.add_option('-n','--normalize',dest='normalize',help='If option is given, all plots will be normalized to the number of entries in the first histogram given in the input list.',action='store_true',default=False)
   parser.add_option('-e','--errors',dest='calc_errors',help='If this option is given, the statistical uncertainties will be shown.',action='store_true',default=False)
   parser.add_option('-y','--logy-pages',dest='logy_pages',help='comma separated list of the page numbers in which to use a log scale on the y-axis.')
   parser.add_option('-x','--logx-pages',dest='logx_pages',help='comma separated list of the page numbers in which to use a log scale on the x-axis.')
   parser.add_option('--minx-pages',dest='minx_pages',help='comma separated list of the page numbers in which to set a minimum scale on the x-axis. Example: "1:-1,4:10" would mean on page one set the x-axis minimum scale to -1, and on page 4 set it to 10.')
   parser.add_option('--maxx-pages',dest='maxx_pages',help='comma separated list of the page numbers in which to set a maximum scale on the x-axis. Example: "1:-1,4:10" would mean on page one set the x-axis maximum scale to -1, and on page 4 set it to 10.')
   parser.add_option('--miny-pages',dest='miny_pages',help='comma separated list of the page numbers in which to set a minimum scale on the y-axis. Example: "1:-1,4:10" would mean on page one set the y-axis minimum scale to -1, and on page 4 set it to 10.')
   parser.add_option('--maxy-pages',dest='maxy_pages',help='comma separated list of the page numbers in which to set a maximum scale on the y-axis. Example: "1:-1,4:10" would mean on page one set the y-axis maximum scale to -1, and on page 4 set it to 10.')
   
   options,args = parser.parse_args()

   if options.input is None:
      parser.error('Must specify -i')

   # parse log axes
   page_settings = {}
   if options.logx_pages is not None:
      for page in options.logx_pages.split(','):
         page = int(page)
         ps = page_settings.get(page)
         if ps is None:
            ps = PageSettings(page)
            page_settings[page] = ps
         ps.log_x = True
   if options.logy_pages is not None:
      for page in options.logy_pages.split(','):
         page = int(page)
         ps = page_settings.get(page)
         if ps is None:
            ps = PageSettings(page)
            page_settings[page] = ps
         ps.log_y = True
   minx_pages = []
   if options.minx_pages is not None:
      for page in options.minx_pages.split(','):
         split = page.split(':')
         page_number = int(split[0])
         value = float(split[1])
         ps = page_settings.get(page_number)
         if ps is None:
            ps = PageSettings(page_number)
            page_settings[page_number] = ps
         ps.min_x = value
   maxx_pages = []
   if options.maxx_pages is not None:
      for page in options.maxx_pages.split(','):
         split = page.split(':')
         page_number = int(split[0])
         value = float(split[1])
         ps = page_settings.get(page_number)
         if ps is None:
            ps = PageSettings(page_number)
            page_settings[page_number] = ps
         ps.max_x = value
   miny_pages = []
   if options.miny_pages is not None:
      for page in options.miny_pages.split(','):
         split = page.split(':')
         page_number = int(split[0])
         value = float(split[1])
         ps = page_settings.get(page_number)
         if ps is None:
            ps = PageSettings(page_number)
            page_settings[page_number] = ps
         ps.min_y = value
   maxy_pages = []
   if options.maxy_pages is not None:
      for page in options.maxy_pages.split(','):
         split = page.split(':')
         page_number = int(split[0])
         value = float(split[1])
         ps = page_settings.get(page_number)
         if ps is None:
            ps = PageSettings(page_number)
            page_settings[page_number] = ps
         ps.max_y = value
   
   # turn off graphics
   ROOT.gROOT.SetBatch(True)
   # take over the ROOT logger
   try:
      root_logger = logging.getLogger('ROOT')
      root_logger.addHandler(logging._handlers.keys()[0])
   except:
      pass
   # load the ATLAS Style rules
   import root.AtlasStyle as AtlasStyle
   AtlasStyle.SetAtlasStyle()
   
   # make file list
   file_list = options.input.split(',')

   # need more than one file to continue
   if len(file_list) < 2:
      raise Exception('Doh! You gotta pass more than one file to make a comparison, right?')

   # loop by files and read in all historams
   rootFiles = []
   for filename in file_list:
      if os.path.exists(filename):
         rf = RootFile(filename)
         rootFiles.append(rf)
      else:
         raise Exception('File does not exist: ' + filename)


   # loop by histogram and compare
   comparison_histos = {}
   for main_hist in rootFiles[0].th1.values():
      main_int = main_hist.Integral('width')
      if options.calc_errors:
         main_hist.Sumw2()
      comparison_histos[main_hist.GetName()] = []
      
      for file in rootFiles[1:]:
         comp_hist = main_hist.Clone(main_hist.GetName() + '_compared')
         comp_hist.SetMarkerStyle(1)
         
         sec_hist = file.th1.get(main_hist.GetName())
         if sec_hist is None:
            logger.warning('Did not find histogram ' + main_hist.GetName() + ' in file ' + file.filename + ' so skipping it.')

         if options.calc_errors:
            sec_hist.Sumw2()

         if options.normalize:
            sec_hist = sec_hist.Clone(sec_hist.GetName() + '_norm')
            scale = main_int/sec_hist.Integral('width')
            sec_hist.Scale(scale)

         comp_hist.Add(sec_hist,-1)
         comp_hist.Divide(sec_hist)

         # extract chi-squared
         chisq = get_chi_square(comp_hist)
         comp_hist.GetYaxis().SetTitle(' (X - Y) / X   (#chi^{2}/ndf = ' + str(chisq) + ')')

         comparison_histos[main_hist.GetName()].append(comp_hist)

   canFile = CanvasFile()

   draw_options = ''
   if options.calc_errors:
      draw_options = 'E1'

   page_number = 0
   for histo_key in sorted(comparison_histos):
      list = comparison_histos[histo_key]
      page_number += 1
      ps = PageSettings(page_number)
      if page_number in page_settings.keys():
         ps = page_settings[page_number]
      canFile.AddPage(log_x = ps.log_x,
                      log_y = ps.log_y,
                      min_x = ps.min_x,
                      max_x = ps.max_x,
                      min_y = ps.min_y,
                      max_y = ps.max_y,
                     )
      for plot in list:
         canFile.AddPlot(plot,draw_options)

   canFile.SaveAs(options.output)



if __name__ == "__main__":
   logging.basicConfig(level=logging.INFO)
   main()
