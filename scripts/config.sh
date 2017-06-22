#!/bin/bash
# Cmake config for biotissiue
# usage ./config.sh [build_type] [logrun_flag]
source $DEVROOT/scripts/util
ROOT=$DEVROOT/bio
LOGRUN_OVERRIDE=$2
if [ -z $1 ]; then
  BUILD_TYPE=Debug
else
  BUILD_TYPE=$1
fi
if [ "$BUILD_TYPE" == "Debug" ] ; then
  BUILD_DIR=$ROOT/build_debug
elif [ "$BUILD_TYPE" == "Release" ] ; then
  BUILD_DIR=$ROOT/build_release
fi
LOGRUN=ON
if [ "$LOGRUN_OVERRIDE" != "" ] ; then
  LOGRUN=$LOGRUN_OVERRIDE
fi
verify_directory_recreate ${BUILD_DIR}
cd $BUILD_DIR
module load cmake
HOSTNAME=`hostname`
if [ "$HOSTNAME" == "q.ccni.rpi.edu" ]; then
  cmake \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
    -DCMAKE_C_COMPILER="mpicc" \
    -DCMAKE_CXX_COMPILER="mpicxx" \
    -DCMAKE_INSTALL_PREFIX=$DEVROOT/install/bio/ \
    -DBUILD_TESTS=OFF \
    -DLOGRUN=$LOGRUN \
    -DSIM_MPI=bgmpi \
    -DCMAKE_PREFIX_PATH=$DEVROOT/install/amsi/lib/cmake/amsi \
    -DSCOREC_DIR=$DEVROOT/install/core/lib/cmake/SCOREC \
    -DSPARSKIT_DIR=$DEVROOT/install/sparskit \
    -DSCORECUTIL_DIR=$DEVROOT/install/scorecutil/xl/ \
    ..
  chmod g+rw $BUILD_DIR
else
  module load $DEVROOT/module/openmpi/1.10.0
  CC=`which mpicc`
  CXX=`which mpicxx`
  cmake \
  -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
  -DBUILD_TESTS=ON \
  -DCMAKE_C_COMPILER=$CC \
  -DCMAKE_CXX_COMPILER=$CXX \
  -DLOGRUN=TRUE \
  -DCMAKE_INSTALL_PREFIX=$DEVROOT/install/bio/ \
  -DCMAKE_PREFIX_PATH=$DEVROOT/install/amsi/lib/cmake/amsi \
  -DSIM_MPI=openmpi110 \
  -DSCOREC_DIR=$DEVROOT/install/core/lib/cmake/SCOREC \
  -DSPARSKIT_DIR=$DEVROOT/install/sparskit/ \
  -DSCORECUTIL_DIR=$DEVROOT/install/scorecutil/openmpi-1.10.0/ \
  ..
fi



