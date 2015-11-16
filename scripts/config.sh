#!/bin/bash

module load simmetrix/simModSuite

HOSTNAME=`hostname`

if [ -z $1 ]; then
  BUILD_TYPE=Debug
else
  BUILD_TYPE=$1
fi

if [ "$BUILD_TYPE" == "Debug" ]; then
  BUILD_DIR=../build_debug
elif [ "$BUILD_TYPE" == "Release" ] ; then
  BUILD_DIR=../build_release
fi

if [ ! -d $BUILD_DIR ]; then
  mkdir $BUILD_DIR
fi
cd $BUILD_DIR
rm -rf ./* #stupid and dangerous

if [ "$HOSTNAME" == "q.ccni.rpi.edu" ]; then
  export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$DEVROOT/install

  cmake \
  -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
  -DBUILD_TESTS=OFF \
  -DLOGRUN=TRUE \
  -DCMAKE_INSTALL_PREFIX=$DEVROOT/install \
  -DSIM_MPI=openmpi14 \
  -DCMAKE_TOOLCHAIN_FILE=$CMAKE_XL_TOOLCHAIN \
  ..
else
  cmake \
  -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
  -DBUILD_TESTS=ON \
  -DLOGRUN=TRUE \
  -DCMAKE_PREFIX_PATH=$DEVROOT/install/amsi/sim/openmpi-1.3.3/lib/cmake/amsi \
  -DCMAKE_INSTALL_PREFIX=$DEVROOT/install/biotissue/git/openmpi-1.3.3 \
  -DCORE_DIR=$DEVROOT/install/core-sim/openmpi-1.3.3/ \
  -DSIMWRAPPER_DIR=$DEVROOT/simPartitionWrapper/PartitionWrapper/lib \
  -DSIM_MPI=openmpi14 \
  ..
fi



