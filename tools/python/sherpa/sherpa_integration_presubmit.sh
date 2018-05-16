#!/usr/bin/env bash

# send output to files
LOG_FILE=/tmp/sherpa_tarball/build.log
exec 1> $LOG_FILE 2>&1
echo $(basename $0) START $SECONDS  $(date)
HPC_PATH=/users/hpcusers
EXE_PATH=$HPC_PATH/argobalsam/production/exe
INSTALL_PATH=$(readlink -f $EXE_PATH/Sherpa-install)
JOB_PATH=$PWD
TARBALL=$1

echo INST: $INSTALL_PATH
echo JOB: $JOB_PATH

CVMFS_GCC463=/cvmfs/atlas.cern.ch/repo/sw/atlas-gcc/463/x86_64/slc6/x86_64-slc6-gcc46-opt

export LD_LIBRARY_PATH=$CVMFS_GCC463/lib64:$HPC_PATH/tools/OpenLoops/lib:/lib:/lib64:/usr/lib:/usr/lib64:/usr/local/lib:/usr/local/lib64:/usr/lib64/openmpi/lib:$INSTALL_PATH/lib
export PATH=$CVMFS_GCC463/bin:/usr/lib64/openmpi/bin:/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:$INSTALL_PATH/bin

tar zxf $TARBALL
RETURN_CODE=$?
if [[ $RETURN_CODE != 0 ]]; then
   echo "tarball extract Failed with exit code: " $RETURN_CODE
   exit $RETURN_CODE
fi
$EXE_PATH/find_and_replace_string.py -i SHERPA_INSTALL_PATH -o $INSTALL_PATH -f makelibs
RETURN_CODE=$?
if [[ $RETURN_CODE != 0 ]]; then
   echo "find_and_replace_string install-path makelibs Failed with exit code: " $RETURN_CODE
   exit $RETURN_CODE
fi
$EXE_PATH/find_and_replace_string.py -i SHERPA_INSTALL_PATH -o $INSTALL_PATH -r -p $JOB_PATH/Process
RETURN_CODE=$?
if [[ $RETURN_CODE != 0 ]]; then
   echo "find_and_replace_string install-path Process Failed with exit code: " $RETURN_CODE
   exit $RETURN_CODE
fi
$EXE_PATH/find_and_replace_string.py -i SHERPA_JOB_PATH -o $JOB_PATH -r -p $JOB_PATH/Process
RETURN_CODE=$?
if [[ $RETURN_CODE != 0 ]]; then
   echo "find_and_replace_string job-path Process Failed with exit code: " $RETURN_CODE
   exit $RETURN_CODE
fi
echo 'makelibs onelib'
sh -c ". ./makelibs -m -j 10"
RETURN_CODE=$?
if [[ $RETURN_CODE != 0 ]]; then
   echo "makelibs Failed with exit code: " $RETURN_CODE
   exit $RETURN_CODE
fi

echo 'makelibs proclibs'
sh -c ". ./makelibs -j 10"
RETURN_CODE=$?
if [[ $RETURN_CODE != 0 ]]; then
   echo "makelibs Failed with exit code: " $RETURN_CODE
   exit $RETURN_CODE
fi

echo $(basename $0) END $SECONDS  $(date)
