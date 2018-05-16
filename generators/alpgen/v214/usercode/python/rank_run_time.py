#!/usr/bin/env python
from ROOT import TH1D,TH2D,gROOT,kBlack
gROOT.SetBatch(True)
from CanvasFile import CanvasFile
import AtlasStyle
AtlasStyle.SetAtlasStyle()
import fnmatch,time,os,sys

base_path = '/users/hpcusers/svn/jobpayloads/'
outdir_filename = 'outdir.txt'
alpgen_wgt_filename_pattern = 'alpout.????????.wgt'
rank_start_index = 7
rank_end_index   = 15
this_year = 2014
event_size_bytes = 57.

time_nbins = 15
time_min   = 15
time_max   = 30
events_nbins = 150
events_min   = 0
events_max   = 15000


def main():
   jobs_to_compare = [
                      job_info('13970511800582820','Apr 11 00:23:46 2014','Apr 11 01:01:23 2014',1024*32,'1024A'),
                      job_info('13970511845948010','Apr 10 23:33:44 2014','Apr 11 00:21:23 2014',1024*32,'1024B'),
                      job_info('13970510806387920','Apr 11 12:48:34 2014','Apr 11 13:43:59 2014',2048*32,'2048A'),
                      job_info('13970510856543650','Apr 11 07:11:47 2014','Apr 11 08:06:38 2014',2048*32,'2048B'),
                     ]

   
   canvasFile = CanvasFile()

   canvasFile.AddPage(auto_max_y=True,
                      legend_x1=0.65,
                      legend_x2=0.95,
                      legend_y1=0.95,
                      legend_y2=0.7)


   canvasFile.AddPage(auto_max_y=True,
                      legend_x1=0.65,
                      legend_x2=0.95,
                      legend_y1=0.95,
                      legend_y2=0.7)

   for i in range(len(jobs_to_compare)):
      job = jobs_to_compare[i]
      filename = os.path.join(base_path,job.job_id,outdir_filename)
      ranks = outdir_parser(filename)
      for rank in ranks:
         job.Fill(rank)
      
      job.h_run_time.SetLineColor(kBlack + i)
      job.h_run_time.SetLineWidth(len(jobs_to_compare)-i)
      canvasFile.AddPlot(job.h_run_time,page_number=0,legend_label=str(job.legend_label),legend_options='l')
      
      job.h_events.SetLineColor(kBlack + i)
      job.h_events.SetLineWidth(len(jobs_to_compare)-i)
      canvasFile.AddPlot(job.h_events,page_number=1,legend_label=str(job.legend_label),legend_options='l')


      canvasFile.AddPage()
      canvasFile.AddPlot(job.h_time_vs_events,'colz')
      canvasFile.AddPage()
      canvasFile.AddPlot(job.h_time_vs_rank,'colz')
      canvasFile.AddPage()
      canvasFile.AddPlot(job.h_events_vs_rank,'colz')
      
      job.ranks = ranks
   
   # plot the difference in rank times:
   h_diff_rank_time = TH1D('diff_rank_time','; Run Time (t) Difference For Fixed Rank M: t_{jobX} - t_{jobY}',10,-5,5)
   h_diff_rank_events_job01 = TH1D('diff_rank_events_job01','; Event Output (N) Difference for Fixed Rank M: N_{1024A} - N_{1024B} ',100,-100,100)
   h_diff_rank_events_job02 = TH1D('diff_rank_events_job02','; Event Output (N) Difference for Fixed Rank M: N_{1024A} - N_{2048A} ',100,-10000,10000)
   for rank in range(1024*32):
      reference = jobs_to_compare[0]
      for i in range(1,len(jobs_to_compare)):
         current = jobs_to_compare[i]
         time_diff = reference.ranks[rank].run_time - current.ranks[rank].run_time
         h_diff_rank_time.Fill(time_diff)
         
      events_diff = jobs_to_compare[0].ranks[rank].numevents - jobs_to_compare[1].ranks[rank].numevents
      h_diff_rank_events_job01.Fill(events_diff)
      events_diff = jobs_to_compare[0].ranks[rank].numevents - jobs_to_compare[2].ranks[rank].numevents
      h_diff_rank_events_job02.Fill(events_diff)

      #print events_diff,reference.ranks[rank].numevents,current.ranks[rank].numevents

   
   canvasFile.AddPage(auto_max_y=True)
   canvasFile.AddPlot(h_diff_rank_time)

   
   canvasFile.AddPage(auto_max_y=True,log_y=True)
   canvasFile.AddPlot(h_diff_rank_events_job01)
    
   canvasFile.AddPage(auto_max_y=True,log_y=True)
   canvasFile.AddPlot(h_diff_rank_events_job02)
   

   h_diff_rank_events_job23 = TH1D('diff_rank_events_job03','; Event Output (N) Difference for Fixed Rank M: N_{2048A} - N_{2048B} ',100,-10000,10000)
   for rank in range(2048*32):
      events_diff = jobs_to_compare[2].ranks[rank].numevents - jobs_to_compare[3].ranks[rank].numevents
      h_diff_rank_events_job23.Fill(events_diff)
   
   canvasFile.AddPage(auto_max_y=True,log_y=True)
   canvasFile.AddPlot(h_diff_rank_events_job23)
      

   canvasFile.SaveAs('test.ps')

   return 0

def outdir_parser(filename):
   file = open(filename)
   
   ranks = []

   for line in file:
      words = line.split()
      try:
         match = fnmatch.fnmatch(words[8],alpgen_wgt_filename_pattern)
      except:
         #print 'Match failed.. continuing'
         continue
      if match:
         filename = words[8]
         rank_number = int(filename[rank_start_index:rank_end_index])
         rank = rank_info(rank_number,words[4],words[5],words[6],words[7])
         ranks.append(rank)

   return ranks
               

class rank_info:
   def __init__(self,rank_number,file_size,month,day,time_of_day):
      self.rank_number  = int(rank_number)
      self.file_size    = int(file_size)
      self.month        = month
      self.day          = day
      self.time         = time_of_day
      self.date_time    = ('%s %s %s %s') % (self.day,self.month,this_year,self.time)
      try:
         struct = time.strptime(self.date_time,'%d %b %Y %H:%M')
      except:
         print 'Failed to parse time',sys.exc_info()[1]
      try:
         self.minutes      = time.mktime(struct)/60. + 5*60 # job times are in UTC, rank time is in CST
      except:
         print 'Failed to convert time to minutes:',sys.exc_info()[1]

      self.numevents    = self.file_size/event_size_bytes

      #print self.rank_number,self.date_time,self.minutes

class job_info:
   def __init__(self,job_id,start_time,end_time,num_ranks,legend_label):
      self.job_id = job_id
      self.start_time = start_time
      self.end_time   = end_time
      self.num_ranks = num_ranks
      self.legend_label = legend_label

      self.start_minutes = time.mktime(time.strptime(self.start_time,'%b %d %H:%M:%S %Y'))/60.
      self.end_minutes   = time.mktime(time.strptime(self.end_time,'%b %d %H:%M:%S %Y'))/60.
      self.duration = self.end_minutes - self.start_minutes

      #print self.job_id,self.start_time,self.start_minutes,self.end_time,self.end_minutes,self.duration

      self.h_run_time = TH1D('run_time_'+self.job_id,';Rank Run Time (minutes)',time_nbins,time_min,time_max)
      self.h_events = TH1D('events_'+self.job_id,';Events generated per rank',events_nbins,events_min,events_max)
      self.h_time_vs_events = TH2D('rt_vs_events_'+self.job_id,';Rank Run Time (minutes);Rank Output Events',time_nbins,time_min,time_max,events_nbins,events_min,events_max)
      self.h_time_vs_rank = TH2D('rt_vs_rank_'+self.job_id,';Rank Run Time (minutes);Rank Number',time_nbins,time_min,time_max,self.num_ranks,0,self.num_ranks-1)
      self.h_events_vs_rank = TH2D('events_vs_rank_'+self.job_id,';Events Generated;Rank Number',events_nbins,events_min,events_max,self.num_ranks,0,self.num_ranks-1)

   def Fill(self,rank):
      rank.run_time = rank.minutes - self.start_minutes
      #print time,rank.minutes,self.start_minutes
      self.h_run_time.Fill(rank.run_time)
      self.h_events.Fill(rank.numevents)
      self.h_time_vs_events.Fill(rank.run_time,rank.numevents)
      self.h_time_vs_rank.Fill(rank.run_time,rank.rank_number)
      self.h_events_vs_rank.Fill(rank.numevents,rank.rank_number)
      #print rank.rank_number


if __name__ == '__main__':
   sys.exit(main())


