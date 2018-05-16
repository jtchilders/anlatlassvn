import ROOT

class CanvasFile:
   
   class PagePlot:
      def __init__(self,plot,draw_options='',legend_label=None,legend_options=None):
         self.plot = plot
         self.draw_options = draw_options
         self.legend_label = legend_label
         self.legend_options = legend_options
      def draw(self,min_x=None,max_x=None,min_y=None,max_y=None,same=False):
         #self.plot.Print("all")
         self.plot_clone = self.plot.Clone(self.plot.GetName()+'_clone')
         if max_y is not None:
            self.plot_clone.SetMaximum(max_y)
         if min_y is not None:
            self.plot_clone.SetMinimum(min_y)
         if min_x is not None:
            self.plot_clone.GetXaxis().SetLimits(min_x,self.plot_clone.GetXaxis().GetXmax())
         if max_x is not None:
            self.plot_clone.GetXaxis().SetLimits(self.plot_clone.GetXaxis().GetXmin(),max_x)
         draw_opts = self.draw_options
         if same:
            draw_opts += 'same'
         #print 'draw opts:',draw_opts
         self.plot_clone.Draw(draw_opts)
   
   class Page:
      def __init__(self,
                   log_x=False,
                   log_y=False,
                   max_y=None,
                   min_y=None,
                   max_x=None,
                   min_x=None,
                   legend_x1 = None, # x is left to right, unit = fraction of current TPad/TCanvas
                   legend_x2 = None, # x1 is the left side, x2 is the right side
                   legend_y1 = None, # y is top to bottom, unit = fraction of current TPad/TCanvas
                   legend_y2 = None, # y1 is the top side, y2 is the bottom side
                   auto_max_y = False,
                   auto_min_y = False,
                   auto_max_x = False,
                   auto_min_x = False,
                  ):
         self.plots = []
         self.log_x = log_x
         self.log_y = log_y
         self.max_y = max_y
         self.min_y = min_y
         self.max_x = max_x
         self.min_x = min_x
         self.legend = None
         if ( legend_x1 is not None and
              legend_x2 is not None and
              legend_y1 is not None and
              legend_y2 is not None):
            self.legend = ROOT.TLegend(legend_x1,legend_y1,legend_x2,legend_y2)
            self.legend.SetFillStyle(0)
            self.legend.SetBorderSize(0)
         self.auto_max_y = auto_max_y
         self.auto_min_y = auto_min_y
         self.auto_max_x = auto_max_x
         self.auto_min_x = auto_min_x
      def add_plot(self,plot,draw_options='',legend_label=None,legend_options=None):
         self.plots.append(CanvasFile.PagePlot(plot,draw_options,legend_label,legend_options))
         if self.legend is not None:
            self.legend.AddEntry(plot,legend_label,legend_options)
      def draw(self):
         max_y = self.max_y
         if self.auto_max_y:
            max_y = get_max_y(self.plots)
         for i in range(len(self.plots)):
            plot = self.plots[i]
            if i == 0:
               plot.draw(self.min_x,self.max_x,self.min_y,max_y)
            else:
               plot.draw(self.min_x,self.max_x,self.min_y,max_y,True)
         if self.legend is not None:
            self.legend.Draw('same')
  
   def __init__(self,
                name = 'canvas_file',
                title = 'CanvasFile',
                position_x = 0, # x is left to right, 0 = left-side, unit = pixels
                position_y = 0, # y is top to bottom, 0 = top-side, unit = pixels
                width = 800, # horizontal width in pixels
                height = 600, # vertical width in pixels
               ):
      self.canvas = ROOT.TCanvas(name,title,position_x,position_y,width,height)
      self.canvas.SetMargin(0.2,0.05,0.15,0.05) # left, right, bottom, top, unit = fraction of pad
      self.pages = []
      self.current_page = -1

   def AddPage(self,
               log_x=False,
               log_y=False,
               max_y=None,
               min_y=None,
               max_x=None,
               min_x=None,
               legend_x1 = None, # x is left to right, unit = fraction of current TPad/TCanvas
               legend_x2 = None, # x1 is the left side, x2 is the right side
               legend_y1 = None, # y is top to bottom, unit = fraction of current TPad/TCanvas
               legend_y2 = None, # y1 is the top side, y2 is the bottom side
               auto_max_y=False,
               auto_min_y=False,
               auto_max_x=False,
               auto_min_x=False,
              ):
      self.pages.append(CanvasFile.Page(log_x,log_y,max_y,min_y,max_x,min_x,
                                        legend_x1,legend_x2,legend_y1,legend_y2,
                                        auto_max_y,auto_min_y,auto_max_x,auto_min_x
                                       )
                        )
      self.current_page = len(self.pages)-1

   def AddPlot(self,plot,draw_options = '',page_number=None,legend_label=None,legend_options=None):
      if page_number is None:
         page_number = self.current_page
      if page_number < len(self.pages):
         self.pages[page_number].add_plot(plot,draw_options,legend_label,legend_options)

   def SaveAs(self,filename):
      # loop over pages, drawing and saving each one
      for i in range(len(self.pages)):
         page = self.pages[i]
         self.canvas.cd()
         page.draw()
         save_filename = filename
         if len(self.pages) > 1:
            if i  == 0: # first page must open canvas file
               save_filename += '('
            elif i == len(self.pages) - 1: # last page must close canvas file
               save_filename += ')'
         self.canvas.SetLogx(page.log_x)
         self.canvas.SetLogy(page.log_y)
         self.canvas.SaveAs(save_filename)


def get_max_y(plots):
   max_y = 0.
   for plot in plots:

      if isinstance(plot.plot,ROOT.TH1):
         #print 'max_y:',max_y
         hist = plot.plot
         for bin_num in range(1,hist.GetNbinsX()):
            y = hist.GetBinContent(bin_num)
            #print 'hist bin:',bin_num,y
            if y > max_y:
               max_y = y
      elif isinstance(plot.plot,ROOT.TGraph):
         graph = plot.plot
         x = None
         y = None
         for point_index in range(0,graph.GetN()):
            graph.GetPoint(point_index,x,y)
            if y > max_y:
               max_y = y
   return max_y*1.1


