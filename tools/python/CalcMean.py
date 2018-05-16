
import math,logging,decimal,json
logger = logging.getLogger(__name__)

class CalcMean:
   def __init__(self):
      self.mean = decimal.Decimal(0)
      self.sigma = decimal.Decimal(0)
      self.n = decimal.Decimal(0)
      self.sum = decimal.Decimal(0)
      self.sum2 = decimal.Decimal(0)

   def add_value(self,value):
      tmpv = decimal.Decimal(str(value))
      self.n += decimal.Decimal(1)
      self.sum += tmpv
      self.sum2 += tmpv*tmpv
      self.mean = decimal.Decimal(0)
      self.sigma = decimal.Decimal(0)

   def calc_mean(self):
      if self.mean != 0:
         return self.mean
      if self.n == 0:
         return decimal.Decimal(0)

      self.mean = self.sum/self.n
      return self.mean

   def calc_sigma(self):
      if self.sigma != 0:
         return self.sigma
      if self.n <= 1:
         return decimal.Decimal(0)
      mean = self.calc_mean()
      try:
         self.sigma = decimal.Decimal(str(math.sqrt( self.sum2/decimal.Decimal(self.n) - mean*mean)))
      except ValueError,e:
         logger.error(str(e) + ': n=' + str(self.n) + ', sum2=' + str(self.sum2) + ', mean=' + str(mean) + ', sum2/n - mean*mean = ' + str(self.sum2/self.n - mean*mean) + ' sum2/n = ' + str(self.sum2/self.n))
         raise
      return self.sigma

   def __add__(self,other):
      new = CalcMean()
      new.n = self.n + other.n
      new.sum = self.sum + other.sum
      new.sum2 = self.sum2 + other.sum2
      return new

   def __eq__(self,other):
      if self.calc_mean() == other.calc_mean():
         return True
      return False

   def __ne__(self,other):
      if self.calc_mean() != other.calc_mean():
         return True
      return False

   def __gt__(self,other):
      if self.calc_mean() > other.calc_mean():
         return True
      return False

   def __lt__(self,other):
      if self.calc_mean() < other.calc_mean():
         return True
      return False

   def __ge__(self,other):
      if self.calc_mean() >= other.calc_mean():
         return True
      return False

   def __le__(self,other):
      if self.calc_mean() <= other.calc_mean():
         return True
      return False


   def get_string(self,format = '%f +/- %f',
                       show_percent_error=False,
                       show_percent_error_format = '%12.2f +/- %12.2f (%5.2f%%)'
                 ):
      s = ''
      if show_percent_error and self.calc_mean() != 0:
         percent_error = self.calc_sigma()/self.calc_mean()*decimal.Decimal(100)
         s = show_percent_error_format % (self.calc_mean(),self.calc_sigma(),percent_error)
      else:
         s = format % (self.calc_mean(),self.calc_sigma())
      return s

class CalcMeanEncoder(json.JSONEncoder):
   def default(self, obj):
      if isinstance(obj, CalcMean):
         return {'mean':float(obj.calc_mean()),'sigma':float(obj.calc_sigma()),'n':int(obj.n)}
      # Let the base class default method raise the TypeError
      return json.JSONEncoder.default(self, obj)
