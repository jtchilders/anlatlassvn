#!/usr/bin/env python
from AlpgenInputFile import AlpgenInputFile as AIF

import logging,glob
logger = logging.getLogger(__name__)


f = AIF()
f.read('template_wenuccjets.txt')

zflave = ['e','mu','tau']

for flav in zflave:
   f.options['iwdecmod'].value = zflave.index(flav) + 1

   for njets in range(0,6):
      f.options['njets'].value = njets
      f.write('w%snucc%ijets.txt' % (flav,njets))



