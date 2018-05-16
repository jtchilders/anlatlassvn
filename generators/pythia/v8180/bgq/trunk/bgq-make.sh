#!/usr/bin/env bash
PATH_TO_HEPMC=/users/hpcusers/svn/tools/HepMC-2.06.08
COPTS_OPTIONAL=
LOPTS_OPTIONAL=
if [[ "$1" == "debug" ]];
then
   echo "BlueGeneQ Make: Setting debug compiler option"
   COPTS_OPTIONAL=-g
fi
if [[ "$1" == "profile" ]];
then
   echo "BlueGeneQ Make: Setting profile compiler option"
   COPTS_OPTIONAL=-pg
   LOPTS_OPTIONAL=-pg
fi
echo $COPTS_OPTIONAL
echo $LOPTS_OPTIONAL
export USRLDFLAGSSHARED=$LOPTS_OPTIONAL
export USRCXXFLAGS=$COPTS_OPTIONAL
./configure --with-hepmc=$PATH_TO_HEPMC --enable-shared 
make

