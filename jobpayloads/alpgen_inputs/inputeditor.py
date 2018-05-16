#!/usr/bin/env python
from AlpgenInputFile import AlpgenInputFile as AIF,AlpgenOption as AO

import logging,glob
logger = logging.getLogger(__name__)

l = glob.glob('input_cards/*.dat')

for file in l:
   f = open(file)
   aif = AIF()
   for line in f:
      try:
         x = int(line[:4])
      except: continue
      
      if 2 <= x and x <= 200 and x != 4:
         parts = line[:-1].split()
         if len(parts) == 4:
            value = float(parts[1])
            name = parts[3]
            
            if value > 0:
               if name == 'izdecmod': name = 'izdecmode'
               aif.options[name] = AO(name,value)
      if x == 0:
         parts = line[:-1].split()
         m = parts[7].split(',')
         for name in m:
            aif.options[name] = AO(name,parts[m.index(name)])
         



   if len(aif.options) > 0:
      aif.write(file.replace('.dat','.txt'))




