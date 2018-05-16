#!/usr/bin/env python
from AlpgenInputFile import AlpgenInputFile as AIF

import logging,glob,os
logger = logging.getLogger(__name__)

l = glob.glob('/users/hpcusers/svn/jobpayloads/alpgen_inputs/vjets_13TeV/*.txt')

for file in l:
   if 'template' in file: continue
   f = AIF()
   f.read(file)
   basefile = os.path.basename(file)
   
   f.options['qfac'].value = 2
   f.write(basefile.replace('.txt','_qfac2.txt'))
   
   f.options['qfac'].value = 0.5
   f.write(basefile.replace('.txt','_qfac0p5.txt'))

   f.options['qfac'].value = 1
   f.options['ktfac'].value = 2
   f.write(basefile.replace('.txt','_ktfac2.txt'))

   f.options['ktfac'].value = 0.5
   f.write(basefile.replace('.txt','_ktfac0p5.txt'))




