###############################################################
#
# Job options file
#
#==============================================================

###############################################################
###############################################################
__doc__ = """
Note: This job options assumes you have run the following commands:
 1) create an EVGEN file (typically done in a previous/next release
athena -c'EVTMAX=10' McParticleTests/iotest_WriteGenEvent_jobOptions.py

 2) create an ASCII-EVGEN file
athena -c'OUTPUT="/tmp/hepmc.1.ascii"' McParticleAlgs/GenEventAsciiWriter_jobOptions.py

 3) create another ASCII-EVGEN file
athena -c'OUTPUT="/tmp/hepmc.2.ascii"' McParticleAlgs/GenEventAsciiWriter_jobOptions.py

 4) read all these ASCII files (ie: run this jobO)
athena McAsciiEventSelector/Example_McAsciiReader_jobOptions.py 
"""
###############################################################
###############################################################

#--------------------------------------------------------------
# General Application Configuration options
#--------------------------------------------------------------
from AthenaCommon.AppMgr import (ServiceMgr as svcMgr,
                                 theApp)

## configure the special Ascii event selector
import McAsciiEventSelector.ReadMcAscii

###################
## User parameters
if not 'MCEVENTKEY' in dir():
    MCEVENTKEY = 'GEN_EVENT'
    pass

if not 'streamName' in dir():
   streamName = 'StreamEVGEN'

if not 'outputFilename' in dir():
   outputFilename = 'myoutput.pool.root'

if not 'ecmEnergy' in dir():
   ecmEnergy = 13000

if not 'runNumber' in dir():
   runNumber = 999999

if not 'inputFilename' in dir():
   # if no input given use this default
   inputFilename = 'myinput.hepmc'
   if 'inputGlob' in dir():
      import glob
      inputCollections = glob.glob(inputGlob)
   # if no input collection given create one
   if not 'inputCollections' in dir():
      inputCollections = [ inputFilename ]

else:
   inputCollections = [ inputFilename ]

#--------------------------------------------------------------
# Event related parameters
#--------------------------------------------------------------
if not 'EVTMAX' in dir():
    EVTMAX = -1

theApp.EvtMax = EVTMAX

svcMgr.EventSelector.InputCollections = inputCollections
svcMgr.EventSelector.OutputLevel = VERBOSE
svcMgr.EventSelector.RunNumber = runNumber
svcMgr.McAsciiCnvSvc.OutputLevel = VERBOSE
svcMgr.McAsciiCnvSvc.McEventsOutput = MCEVENTKEY

import EventInfoMgt.EventInfoMgtInit

svcMgr.TagInfoMgr.ExtraTagValuePairs += ["beam_energy", str(int(ecmEnergy*Units.GeV/2.0))]
svcMgr.TagInfoMgr.ExtraTagValuePairs += ["beam_type", 'collisions']

## Add special config option (extended model info for BSM scenarios)
svcMgr.TagInfoMgr.ExtraTagValuePairs += ["specialConfiguration", 'NONE' ]

## update atlas release
import os
project = os.environ ['AtlasProject']
version = os.environ ['AtlasPatchVersion']
release = project + '-' + version
if 'AtlasRelease' in svcMgr.TagInfoMgr.ExtraTagValuePairs:
   index = svcMgr.TagInfoMgr.ExtraTagValuePairs.index('AtlasRelease')
   svcMgr.TagInfoMgr.ExtraTagValuePairs[index+1] = release
else:
   svcMgr.TagInfoMgr.ExtraTagValuePairs += ["AtlasRelease",release]


#--------------------------------------------------------------
# Private Application Configuration options
#--------------------------------------------------------------
# Load "user algorithm"
#top algorithms to be run, and the libraries that house them

from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

########
# Read the McEventCollection
#
if not 'DUMP' in dir():
    DUMP = True
    pass
if DUMP:
    job += CfgMgr.DumpMC( "ReadGenEvent",
                          McEventKey = MCEVENTKEY )

#---------------------------------------------------------------
# Pool Persistency
#---------------------------------------------------------------
from AthenaPoolCnvSvc.WriteAthenaPool import AthenaPoolOutputStream
outStream = AthenaPoolOutputStream(streamName,outputFilename)
outStream.ItemList  = [
    "EventInfo#*",
    "McEventCollection#*",
    ]

#--------------------------------------------------------------
# Set output level threshold (2=DEBUG, 3=INFO, 4=WARNING, 5=ERROR, 6=FATAL )
#--------------------------------------------------------------
#svcMgr.MessageSvc.OutputLevel = DEBUG

#==============================================================
#
# End of job options file
#
###############################################################
