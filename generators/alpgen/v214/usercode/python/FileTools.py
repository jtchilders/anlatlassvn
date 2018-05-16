import subprocess,logging
logger = logging.getLogger(__name__)

def GetLineCount(filename):
   logging.debug('Getting line count for ' + filename)
   p = subprocess.Popen(['wc','-l',filename],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   stdout,stderr = p.communicate()
   return int(stdout.split()[0])




