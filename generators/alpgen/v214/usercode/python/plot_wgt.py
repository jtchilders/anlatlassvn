#!/usr/bin/env python
import logging,os,sys,Queue,math,array,multiprocessing,time,optparse,glob
logger = logging.getLogger(__name__)

from AlpgenWgtFile import AlpgenWgtFile
from CanvasFile import CanvasFile


STATUS_INTERVAL = 1000000
NUM_PROCESSES = 12

def get_log_bins(min,max,nbins):
   edges = []

   if min <= 0 or max <= 0 or max < min:
      logger.error(' bin min and max must be greater than zero: ' + min + ' ' + max )
      sys.exit(-1)

   total_width = math.log(max) - math.log(min)
   #logger.info('width ' + str(total_width))
   if(nbins > 0):
      log_bin_size = total_width/nbins
      nominal_bin_size = math.exp(log_bin_size) - 1
   #logger.info('bin size ' + str(log_bin_size))
   current_x = min
   edges.append(current_x)
   for bin in range(nbins):
      current_x += current_x*nominal_bin_size
      #logger.info(' > '+str(current_x)+' '+str(linear_width))
      edges.append(current_x)
   
   return edges
#get_log_bins(1,10**6,100)
#sys.exit(-1)


def plot_events((filename,histograms,thread_id)):
   try:
      logger = logging.getLogger(__name__+thread_id)
      logger.info('Plotting events in file '+filename)
      wgt_file = AlpgenWgtFile(filename)
      histograms_local = []
      for histogram in histograms:
         histograms_local.append(histogram.Clone(histogram.GetName() + '_' + thread_id))
      event_counter = 0
      total_events = wgt_file.nevts
      one_percent = int(total_events*0.01) + 1
      while wgt_file.LoadEvent(): # and event_counter < max_events:
         if event_counter % one_percent == 0 and total_events > 0:
            logger.info('Percent Done: %3d%%' % int((event_counter + 1 )*100 / total_events))
         histograms_local[0].Fill(wgt_file.iseed1)
         histograms_local[1].Fill(wgt_file.iseed2)
         histograms_local[2].Fill(wgt_file.iavg_store)
         histograms_local[3].Fill(wgt_file.weight)
         if wgt_file.weight > 0:
            histograms_local[4].Fill(1/wgt_file.weight)
         histograms_local[5].Fill(wgt_file.x1)
         event_counter += 1
      wgt_file.close()
      #print '<',thread_id,'> queuing histograms'
      #q.put([h_iseed1,h_iseed2,h_iavg_store,h_weight,h_x1])
      logger.info(' done')

      return histograms_local
   except:
      logger.exception('Exception received.')
      return [None]*len(histograms)



def main():
   logging.basicConfig(level=logging.INFO)
   parser = optparse.OptionParser(description='Plot the data contained in an AlpGen Weighted Event file (*.wgt/*.par).')
   parser.add_option('-i','--input-file',dest='input_file',help='Input filename.')
   parser.add_option('-g','--input-glob',dest='input_glob',help='Input glob pattern to get file list, for example: "/path/file*.wgt". Must use quotations around patter.')
   options,args = parser.parse_args()

   file_list = []
   if options.input_file is not None and os.path.exists(options.input_file):
      file_list.append(options.input_file)
   if options.input_glob is not None:
      file_list += glob.glob(options.input_glob)

   if len(file_list) == 0:
      parser.error('no input files to process')

   logger.info('looping over ' + str(len(file_list)) + ' files.')
   
   try:
      import ROOT
   except:
      print 'ROOT not in PYTHONPATH'
      sys.exit(-1)
   ROOT.gROOT.SetBatch(True)
   try:
      root_logger = logging.getLogger('ROOT')
      root_logger.addHandler(logging._handlers.keys()[0])
   except:
      pass
   import AtlasStyle
   AtlasStyle.SetAtlasStyle()

   # create histograms
   h_iseed1 = ROOT.TH1D('h_iseed1',';iseed1',100,0,2**31)
   h_iseed2 = ROOT.TH1D('h_iseed2',';iseed2',100,0,2**31)
   h_iavg_store = ROOT.TH1D('h_iavg_store',';iavg_store',25,0,25)
   log_bins = get_log_bins(10**-11,10**9,100)
   h_weight = ROOT.TH1D('h_weight',';weight',len(log_bins)-1,array.array('d',log_bins))
   log_bins = get_log_bins(1e-9,10**11,100)
   h_weight_inverse = ROOT.TH1D('h_weight_inverse',';weight^{-1}',len(log_bins)-1,array.array('d',log_bins))
   log_bins = get_log_bins(1e-4,1,100)
   h_x1 = ROOT.TH1D('h_x1',';x1',len(log_bins)-1,array.array('d',log_bins))

   histograms = [h_iseed1,h_iseed2,h_iavg_store,h_weight,h_weight_inverse,h_x1]

   args = []
   i = 0
   for filename in file_list:
      args.append([filename,histograms,str(i)])
      i += 1
   
   p = multiprocessing.Pool(processes=NUM_PROCESSES)
   logging.info('launch processes')
   #logger.info(str(args))
   try:
      results = p.map(plot_events,args)
   except:
      logger.exception('Exception received in Process Pool')
      raise


   # combine all the histograms
   for histGroup in results:
      if histGroup[0] is None:
         continue
      h_iseed1.Add(histGroup[0])
      h_iseed2.Add(histGroup[1])
      h_iavg_store.Add(histGroup[2])
      h_weight.Add(histGroup[3])
      h_weight_inverse.Add(histGroup[4])
      h_x1.Add(histGroup[5])
      

   can = CanvasFile()
   can.AddPage()
   can.AddPlot(h_iseed1)
   can.AddPage()
   can.AddPlot(h_iseed2)
   can.AddPage(log_y=True)
   can.AddPlot(h_iavg_store)
   can.AddPage(log_x=True,log_y=True)
   can.AddPlot(h_weight)
   can.AddPage(log_x=True,log_y=True)
   can.AddPlot(h_weight_inverse)
   can.AddPage(log_y=True)
   can.AddPlot(h_x1)

   can.SaveAs("plot_wgt.py.ps")
   
   f = ROOT.TFile('plot_wgt.py.root','RECREATE')
   h_iseed1.Write()
   h_iseed2.Write()
   h_iavg_store.Write()
   h_weight.Write()
   h_weight_inverse.Write()
   h_x1.Write()
   f.Close()

   return 0

      




if __name__ == "__main__":
   sys.exit(main())





