#!/usr/bin/env bash
module unload altd
module unload darshan
module unload PrgEnv-intel
module load PrgEnv-gnu
module unload cray-shmem
module load scons
module load mercurial
module load cmake
module load boost
module load gsl

COPTS_OPTIONAL=
LOPTS_OPTIONAL=
INSTALL_PATH=$PWD/local
echo $COPTS_OPTIONAL
echo $LOPTS_OPTIONAL
mkdir -p $INSTALL_PATH
./bootstrap
./configure --with-pic --prefix=$INSTALL_PATH --with-momentum=MEV --with-length=MM CC=cc CXX=CC CFLAGS="$COPPS_OPTIONAL" CXXFLAGS="$COPTS_OPTIONAL" LDFLAGS="$LOPTS_OPTIONAL" --enable-shared=no
make
make install



