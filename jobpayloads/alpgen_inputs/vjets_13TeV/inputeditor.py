#!/usr/bin/env python
from AlpgenInputFile import AlpgenInputFile as AIF

import logging,glob
logger = logging.getLogger(__name__)

l = glob.glob('w*bb*.txt')

for file in l:
   f = AIF()
   f.read(file)
   f.options['ptbmin'].value = 0.0
   f.write(file)



