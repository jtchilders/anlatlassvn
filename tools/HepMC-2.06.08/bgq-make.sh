#!/usr/bin/env bash
#set -o xtrace
COPTS_OPTIONAL=-fPIC
LOPTS_OPTIONAL=
INSTALL_PATH=$PWD/local
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
mkdir -p $INSTALL_PATH
./bootstrap
./configure --prefix=$INSTALL_PATH --with-momentum=MEV --with-length=MM CFLAGS="$COPPS_OPTIONAL" CXXFLAGS="$COPTS_OPTIONAL" LDFLAGS="$LOPTS_OPTIONAL" --with-GENSER=/users/hpcusers/svn/generators/pythia/v6428/unedited --enable-shared=no
make
make install
