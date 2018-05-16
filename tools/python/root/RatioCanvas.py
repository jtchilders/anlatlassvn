from ROOT import TCanvas,TPad,TLegend,TGaxis
import ChiSquared

class Plot:
   def __init__(self,plot,
                draw_options='',
                legend_label=None,
                legend_draw_options='lp',
               ):
      self.plot                  = plot
      self.draw_options          = draw_options
      if legend_label is None:
         self.legend_label       = plot.GetName()
      else:
         self.legend_label       = legend_label
      self.legend_draw_options   = legend_draw_options
   def Draw(self,draw_options = None,same=False):
      if draw_options is not None:
         self.draw_options = draw_options
      if same:
         self.plot.Draw(self.draw_options + ' same')
      else:
         self.plot.Draw(self.draw_options)
   def AddToLegend(self,legend,chi2=None):
      legend_label = self.legend_label
      if chi2 is not None:
         legend_label += ' chi2 = %05.3f' % (chi2/self.plot.GetNbinsX())
      legend.AddEntry(self.plot,legend_label,self.legend_draw_options)
   def CreateRatio(self,reference,ratio_order):
      # create new ratio histogram
      self.ratio = None
      if ratio_order is RatioPad.REFERENCE_IS_DENOMINATOR:
         clone_name = 'ratio_' + self.plot.GetName() + '_' + reference.plot.GetName()
         self.ratio = self.plot.Clone(clone_name)
         self.ratio.Divide(reference.plot)
      else:
         clone_name = 'ratio_' + reference.plot.GetName() + '_' + self.plot.GetName()
         self.ratio = reference.plot.Clone(clone_name)
         self.ratio.Divide(self.plot)
   def DrawRatio(self,draw_options = None, same=False):
      if draw_options is not None:
         self.draw_options = draw_options
      if same:
         self.ratio.Draw(self.draw_options + ' same')
      else:
         self.ratio.Draw(self.draw_options)



   
class TopPad:
   def __init__(self):
      # list of plots for this pad
      self.plots = []

   def AddPlot(self,plot,draw_options='',legend_label=None,legend_draw_options='lp'):
      clone = plot.Clone(plot.GetName() + '_clone')
      clone.SetDirectory(0)
      self.plots.append(Plot(clone,draw_options,legend_label,legend_draw_options))
   
   def Draw(self,
            pad_name             = 'top_pad',
            pad_title            = 'Top Pad',
            pad_x1               = 0,
            pad_y1               = 0,
            pad_x2               = 1,
            pad_y2               = 1,
            pad_border_size      = 1,
            pad_bottom_margin    = 0,
            pad_top_margin       = 0,
            pad_left_margin      = 0,
            pad_right_margin     = 0,
            pad_y_ndiv           = None,
            legend_x1            = None,
            legend_y1            = None,
            legend_x2            = None,
            legend_y2            = None,
            legend_border_size   = 0,
            legend_fill_style    = 0,
            legend_chi2          = False,
            ):

      # create the pad
      self.pad = TPad(pad_name,
                      pad_title,
                      pad_x1,
                      pad_y1,
                      pad_x2,
                      pad_y2,
                      pad_border_size,
                      0,
                     )
      self.pad.SetMargin(pad_left_margin,pad_right_margin,pad_bottom_margin,pad_top_margin)
      self.pad.SetFillColor(0)
      self.pad.Draw()
      self.pad.cd()

      # create the legend
      draw_legend = False
      if legend_x1 != None and legend_x2 != None and legend_y1 != None and legend_y2 != None:
         draw_legend = True

      if draw_legend:
         self.legend = TLegend(legend_x1,legend_y1,legend_x2,legend_y2)
         self.legend.SetBorderSize(legend_border_size)
         self.legend.SetFillStyle(legend_fill_style)

      # draw the plots
      for i in range(len(self.plots)):
         plot = self.plots[i]
         chi2 = None
         if i == 0:
            #plot.plot.GetYaxis().SetLabelSize(0)
            # can turn off the X-axis labels
            plot.plot.GetXaxis().SetLabelSize(0)
            plot.plot.GetXaxis().SetTitle('')
            if pad_y_ndiv is not None:
               plot.plot.GetYaxis().SetNdivisions(pad_y_ndiv)
            plot.Draw()
         else:
            plot.Draw(same=True)
            # calculate chi2 between this plot and the first
            if legend_chi2:
               chi2 = ChiSquared.TH1ChiSquared(self.plots[0].plot,plot.plot)
         if draw_legend:
            plot.AddToLegend(self.legend,chi2)
      
      if draw_legend:
         self.legend.Draw('same')

      # draw y-axis again
      #x_offset = 0
      #histo = self.plots[1].plot
      #y_min = histo.GetMinimum()
      #y_max = histo.GetMaximum()
      #axis = TGaxis(x_offset,y_min,x_offset,y_max,y_min,y_max)
      #axis.SetLabelFont(histo.GetYaxis().GetLabelFont())
      #axis.SetLabelSize(histo.GetYaxis().GetLabelSize())
      #axis.Draw()


class RatioPad:
   REFERENCE_IS_DENOMINATOR = 0
   REFERENCE_IS_NUMERATOR = 1
   def __init__(self): 
      # reference plot: all other plots will be plotted as a ratio to this plot
      self.reference_plot = None

      # list of other plots, which will be calculated as a ratio with the reference
      self.plots = []

   def AddPlot(self,plot,draw_options='',legend_label=None,
                    legend_draw_options='lp',is_reference=False):
      clone = plot.Clone(plot.GetName() + '_clone')
      clone.SetDirectory(0)
      if is_reference:
         self.reference_plot = Plot(clone,draw_options,legend_label,legend_draw_options)
      else:
         self.plots.append(Plot(clone,draw_options,legend_label,legend_draw_options))

   def Draw(self,
            pad_name             = 'bottom_pad',
            pad_title            = 'Bottom Pad',
            pad_x1               = 0,
            pad_y1               = 0,
            pad_x2               = 1,
            pad_y2               = 1,
            pad_border_size      = 1,
            pad_ratio_order      = REFERENCE_IS_DENOMINATOR,
            pad_bottom_margin    = 0,
            pad_top_margin       = 0,
            pad_left_margin      = 0,
            pad_right_margin     = 0,
            yaxis_min            = 0.9,
            yaxis_max            = 1.1,
            y_ndiv               = None,
            y_title              = None,
            ratio_factor         = 2,
           ):
      # create the pad
      self.pad = TPad(pad_name,
                      pad_title,
                      pad_x1,
                      pad_y1,
                      pad_x2,
                      pad_y2,
                      pad_border_size,
                      0
                     )
      self.pad.SetMargin(pad_left_margin,pad_right_margin,pad_bottom_margin,pad_top_margin)
      self.pad.SetFillColor(0)
      self.pad.Draw()
      self.pad.cd()
      
      # create ratio plot for reference
      self.reference_plot.CreateRatio(self.reference_plot,pad_ratio_order)
      # the ratio plot will be plotted first so it determines the axes
      self.reference_plot.ratio.SetMaximum(yaxis_max)
      self.reference_plot.ratio.SetMinimum(yaxis_min)
      # here I setup the ratio to show up 
      self.reference_plot.ratio.GetXaxis().SetLabelSize(self.reference_plot.ratio.GetXaxis().GetLabelSize()*ratio_factor)
      self.reference_plot.ratio.GetYaxis().SetLabelSize(self.reference_plot.ratio.GetYaxis().GetLabelSize()*ratio_factor)
      self.reference_plot.ratio.GetXaxis().SetTitleSize(self.reference_plot.ratio.GetXaxis().GetTitleSize()*ratio_factor)
      self.reference_plot.ratio.GetYaxis().SetTitleSize(self.reference_plot.ratio.GetYaxis().GetTitleSize()*ratio_factor)

      if y_ndiv is not None:
         self.reference_plot.ratio.GetYaxis().SetNdivisions(y_ndiv)
      if y_title is not None:
         self.reference_plot.ratio.GetYaxis().SetTitle(y_title)

      self.reference_plot.DrawRatio()

      for plot in self.plots:
         plot.CreateRatio(self.reference_plot,pad_ratio_order)
         plot.DrawRatio(same=True)
      

      
         
               


class RatioCanvas:
   def __init__(self,
                canvas_name            = 'ratio_canvas',
                canvas_title           = 'Ratio Canvas',
                canvas_x               = 0,
                canvas_y               = 0,
                canvas_width           = 800,
                canvas_height          = 600,
                topPad_name            = 'top_pad',
                topPad_title           = 'Top Pad',
                topPad_x1              = 0, # x coordinate for the bottom-left corner (unit = fraction of Canvas)
                topPad_y1              = 0.3, # y coordinate for the bottom-left corner
                topPad_x2              = 1, # x coord. for the top-right corner
                topPad_y2              = 1, # y coord. for the top-right corner
                topPad_borderSize      = 1,
                topPad_bottomMargin    = 0.03,
                topPad_topMargin       = 0.03,
                topPad_leftMargin      = 0.12,
                topPad_rightMargin     = 0,
                topPad_y_ndiv          = None,
                ratioPad_name          = 'bottom_pad',
                ratioPad_title         = 'Bottom Pad',
                ratioPad_x1            = 0, # bottom-left
                ratioPad_y1            = 0, # bottom-left
                ratioPad_x2            = 1, # top-right
                ratioPad_y2            = 0.3, # top-right
                ratioPad_borderSize    = 1,
                ratioPad_ratioOrder    = RatioPad.REFERENCE_IS_DENOMINATOR,
                ratioPad_bottomMargin  = 0.25,
                ratioPad_topMargin     = 0.03,
                ratioPad_leftMargin    = 0.12,
                ratioPad_rightMargin   = 0,
                ratioPad_yaxis_min     = 0.9,
                ratioPad_yaxis_max     = 1.1,
                ratioPad_y_ndiv        = None,
                ratioPad_y_title       = None,
                legend_x1              = 0.6, # top-left
                legend_y1              = 0.95, # top-left
                legend_x2              = 0.85, # bottom-right
                legend_y2              = 0.7, # bottom-right
                legend_border_size     = 0,
                legend_fill_style      = 0,
                legend_chi2            = False,
                canvas                 = None,
               ):
      # create TopPad and RatioPad
      self.top_pad = TopPad()
      self.ratio_pad = RatioPad()
      self.canvas = canvas

      # Canvas Settings
      self.canvas_name           = canvas_name
      self.canvas_title          = canvas_title
      self.canvas_x              = canvas_x
      self.canvas_y              = canvas_y
      self.canvas_width          = canvas_width
      self.canvas_height         = canvas_height

      # Top Pad settings
      self.topPad_name           = topPad_name
      self.topPad_title          = topPad_title
      self.topPad_x1             = topPad_x1
      self.topPad_y1             = topPad_y1
      self.topPad_x2             = topPad_x2
      self.topPad_y2             = topPad_y2
      self.topPad_borderSize     = topPad_borderSize
      self.topPad_bottomMargin   = topPad_bottomMargin
      self.topPad_topMargin      = topPad_topMargin
      self.topPad_leftMargin     = topPad_leftMargin
      self.topPad_rightMargin    = topPad_rightMargin
      self.topPad_y_ndiv         = topPad_y_ndiv

      # Ratio Pad settings
      self.ratioPad_name         = ratioPad_name
      self.ratioPad_title        = ratioPad_title
      self.ratioPad_x1           = ratioPad_x1
      self.ratioPad_y1           = ratioPad_y1
      self.ratioPad_x2           = ratioPad_x2
      self.ratioPad_y2           = ratioPad_y2
      self.ratioPad_borderSize   = ratioPad_borderSize
      self.ratioPad_ratioOrder   = ratioPad_ratioOrder
      self.ratioPad_bottomMargin = ratioPad_bottomMargin
      self.ratioPad_topMargin    = ratioPad_topMargin
      self.ratioPad_leftMargin   = ratioPad_leftMargin
      self.ratioPad_rightMargin  = ratioPad_rightMargin
      self.ratioPad_ymin         = ratioPad_yaxis_min
      self.ratioPad_ymax         = ratioPad_yaxis_max
      self.ratioPad_y_ndiv       = ratioPad_y_ndiv
      self.ratioPad_y_title      = ratioPad_y_title

      # legend settings
      self.legend_x1             = legend_x1
      self.legend_y1             = legend_y1
      self.legend_x2             = legend_x2
      self.legend_y2             = legend_y2
      self.legend_border_size    = legend_border_size
      self.legend_fill_style     = legend_fill_style
      self.legend_chi2           = legend_chi2

 
   def AddPlot(self,plot,draw_options='',legend_label=None,
               legend_draw_options='lp',is_reference=False):
      self.top_pad.AddPlot(plot,draw_options,legend_label,legend_draw_options)
      self.ratio_pad.AddPlot(plot,draw_options,legend_label,legend_draw_options,is_reference)
  
   def Draw(self,filename=None):
      # create the global canvas
      if( self.canvas is None):
         print 'creating canvas'
         self.canvas = TCanvas(self.canvas_name,
                               self.canvas_title,
                               self.canvas_x,
                               self.canvas_y,
                               self.canvas_width,
                               self.canvas_height,
                              )
         self.canvas.SetMargin(0.,0.,0.,0.)
         self.canvas.SetFillColor(0)
         
      self.canvas.cd()
      # Draw the top pad
      self.top_pad.Draw(
                         self.topPad_name,
                         self.topPad_title,
                         self.topPad_x1,
                         self.topPad_y1,
                         self.topPad_x2,
                         self.topPad_y2,
                         self.topPad_borderSize,
                         self.topPad_bottomMargin,
                         self.topPad_topMargin,
                         self.topPad_leftMargin,
                         self.topPad_rightMargin,
                         self.topPad_y_ndiv,
                         self.legend_x1,
                         self.legend_y1,
                         self.legend_x2,
                         self.legend_y2,
                         self.legend_border_size,
                         self.legend_fill_style,
                         self.legend_chi2,
                       )
      # bring focus back to canvas
      self.canvas.cd()
      # calculate ratio of y-size of top pad and ratio pad
      ratio_factor = (self.topPad_y2 - self.topPad_y1)/(self.ratioPad_y2 - self.ratioPad_y1)
      # create the bottom pad
      self.ratio_pad.Draw(
                          self.ratioPad_name,
                          self.ratioPad_title,
                          self.ratioPad_x1,
                          self.ratioPad_y1,
                          self.ratioPad_x2,
                          self.ratioPad_y2,
                          self.ratioPad_borderSize,
                          self.ratioPad_ratioOrder,
                          self.ratioPad_bottomMargin,
                          self.ratioPad_topMargin,
                          self.ratioPad_leftMargin,
                          self.ratioPad_rightMargin,
                          self.ratioPad_ymin,
                          self.ratioPad_ymax,
                          self.ratioPad_y_ndiv,
                          self.ratioPad_y_title,
                          ratio_factor,
                         )
      if filename is not None:
         self.canvas.SaveAs(filename)
