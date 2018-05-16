#!/usr/bin/env python
from AlpgenInputFile import AlpgenInputFile as AIF

import logging,glob
logger = logging.getLogger(__name__)


f = AIF()
f.read('template_zeejets.txt')

zflave = ['ee','mumu','tautau']

for flav in zflave:
   f.options['izdecmode'].value = zflave.index(flav) + 1

   for njets in range(0,6):
      f.options['njets'].value = njets
      f.write('z%s%ijets.txt' % (flav,njets))



