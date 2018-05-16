#!/usr/bin/env python
from AlpgenInputFile import AlpgenInputFile as AIF

import logging,glob,os,sys
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

filea = sys.argv[1]
fileb = sys.argv[2]

fA = AIF()
fA.read(filea)

fB = AIF()
fB.read(fileb)

diffopts = 0
if len(fA.options) != len(fB.options):
   optsA = set(fA.options.keys())
   optsB = set(fB.options.keys())
   AnotB = optsA - optsB
   BnotA = optsB - optsA
   print 'in A but not B',AnotB
   print 'in B but not A',BnotA


for optname,optvalue in fA.options.iteritems():
   if optname in fB.options:
      if float(optvalue.value) != float(fB.options[optname].value):
         print optname,'differs','A=',optvalue.value,'B=',fB.options[optname].value

