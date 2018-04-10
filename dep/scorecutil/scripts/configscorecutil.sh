#!/bin/bash
# CMake build for SCORECUTIL
#
# usage: ./configscorecutil.sh [build_type] [install_prefix]
#
# todo: check $2 to verify it as a valid install prefix
#

HOSTNAME=`hostname`

if [ "$1" == "Debug" -o "$1" == "Release" ]; then

  cd ../build 

  cmake \
    -DCMAKE_C_COMPILER=mpicc\
    -DCMAKE_CXX_COMPILER=mpicxx\
    -DCMAKE_BUILD_TYPE=$1 \
    -DCMAKE_INSTALL_PREFIX=$2 \
    -DBOOST_INCLUDE_DIR=$DEVROOT/boost_1_56_0/ \
    -DSIM_INCLUDE_DIR=$SIM_INCLUDE_DIR \
    -DSIM_LIB_DIR=$SIM_LIB_DIR \
    ..

fi

make -C ../build clean




