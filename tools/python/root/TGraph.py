import ROOT

class TGraph:
   def __init__(self,name=None):
      self.graph = ROOT.TGraph()
      if name: self.graph.SetName(name)

   def AddPoint(self,x,y):
      self.graph.SetPoint(self.graph.GetN(),x,y)

   def SetAxisLabels(self,x_label=None,y_label=None):
      if x_label:
         self.graph.GetXaxis().SetTitle(x_label)
      if y_label:
         self.graph.GetYaxis().SetTitle(y_label)
