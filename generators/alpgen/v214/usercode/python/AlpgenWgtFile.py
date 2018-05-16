
import logging,sys,os
logger = logging.getLogger(__name__)

class AlpgenWgtFile:
   BYTES_PER_EVENT = 57
   def __init__(self,filename):
      self.filename = filename
      self.file = None
      self.open()
      self.nevts = os.path.getsize(self.filename)/AlpgenWgtFile.BYTES_PER_EVENT

   def open(self):
      try:
         self.file = open(self.filename)
      except:
         logger.error('Cannot Open file '+self.filename+' Exception: '+sys.exec_info()[1])
         return

   def LoadEvent(self,event_num=None):
      if event_num is not None:
         self.file.seek(event_num*AlpgenWgtFile.BYTES_PER_EVENT,0)
      self.line = self.file.readline()
      if len(self.line) == 0:
         return False

      item = self.line.split()
      self.iseed1       = int(  item[0])
      self.iseed2       = int(  item[1])
      self.iavg_store   = float(item[2])
      self.weight       = float(item[3])
      self.x1           = float(item[4])
      return True
   
   def close(self):
      self.file.close()
      self.file = None
