#!/bin/bash
# CMake build for MICRO_FM

module load boost
module load simmetrix/simModSuite
module load apf_1.0
module load FEA_0.1
module load SCORECUTIL_1.0
module load petsc/default

cd ../build 

if [ "$1" == "Debug" -o "$1" == "Release" ]; then

cmake \
  -DCMAKE_BUILD_TYPE=$1 \
  -DCMAKE_INSTALL_PREFIX=$2\
  -DPETSC_DIR=$PETSC_DIR \
  -DBOOST_INCLUDE_DIR=$BOOST_INCLUDE_DIR \
  -DSIM_INCLUDE_DIR=$SIM_INCLUDE_DIR\
  -DSCORECUTIL_INCLUDE_DIR=$SCORECUTIL_INCLUDE_DIR \
  -DFEA_INCLUDE_DIR=$FEA_INCLUDE_DIR \
  -DAPF_INCLUDE_DIR=$APF_INCLUDE_DIR \
  ..

fi

make clean

