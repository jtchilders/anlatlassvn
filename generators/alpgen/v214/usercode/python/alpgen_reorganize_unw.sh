#!/usr/bin/env bash

ARGS=$*
echo ARGS=$ARGS

echo Setup Atlas environment
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh --

echo Setup Rucio
export RUCIO_ACCOUNT=pilot
export X509_USER_PROXY=/users/hpcusers/argobalsam/production/upload_proxy
lsetup rucio

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/users/hpcusers/svn/tools/lhapdf-5.9.1/local/lib

echo Run reorg script
python /users/hpcusers/svn/generators/alpgen/v214/usercode/python/alpgen_reorganize_unw.py $*

