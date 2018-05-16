from ROOT import TH1,TMath
import logging
logger = logging.getLogger(__name__)

def TH1ChiSquared(h1,h2):
   if h1.GetNbinsX() != h2.GetNbinsX():
      logger.exception('TH1ChiSquare: bin count is not equal, h1 = ' + h1.GetNbinsX() + ', h2 = ' + h2.GetNbinsX() )

   chi2 = 0
   bin = 1
   #logger.info(' >> ' + str(h1.GetNbinsX()) + ' ' + str(h2.GetNbinsX()))
   while bin <= h1.GetNbinsX():
      #logger.info('   >> ' + str(bin) + ' ' + str(h1.GetBinContent(bin)) + ' ' + str(h1.GetBinError(bin)) + ' ' + str(h2.GetBinContent(bin)) + ' ' + str(h2.GetBinError(bin)))
      chi2 += chisquared(h1.GetBinContent(bin),h1.GetBinError(bin),h2.GetBinContent(bin),h2.GetBinError(bin))
      bin += 1
   return chi2

def TH1pvalueAndChi2(h1,h2,ndf=None):
   try:
      chi2 = TH1ChiSquare(h1,h2)
   except:
      raise
   pvalue = 0
   if ndf is not None:
      pvalue = TMath.Prob(chi2,ndf)
   else:
      pvalue = TMath.Prob(chi2,h1.GetNbinsX())

   return pvalue,chi2


   
   

# using:
# chi^2 = (x-y)^2 / (dx^2 + dy^2)
def chisquared(x,dx,y,dy):
   try:
      return (x-y)**2 / ( dx**2 + dy**2)
   except ZeroDivisionError:
      return 0

