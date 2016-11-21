#!/bin/bash
# CMake config for MICRO_FO
#
# usage: ./configmicro_fo.sh [build_type] [install_prefix]
#
# todo: check $2 to verify install prefix is valid
#

HOSTNAME=`hostname`

module load Sparskit_1.0
module load simmetrix/simModSuite
module load apf_1.0
module load SCORECUTIL_1.0
module load AMSI_0.1

if [ "$1" == "Debug" -o "$1" == "Release" ]; then

cd ../build

  cmake \
    -DHOSTNAME=$HOSTNAME \
    -DCMAKE_BUILD_TYPE=$1 \
    -DCMAKE_INSTALL_PREFIX=$2 \
    -DBOOST_INCLUDE_DIR=$BOOST_INCLUDE_DIR \
    -DSPARSKIT_INCLUDE_DIR=$SPARSKIT_INCLUDE_DIR \
    -DSIM_INCLUDE_DIR=$SIM_INCLUDE_DIR \
    -DSCORECUTIL_INCLUDE_DIR=$SCOREC_INCLUDE_DIR \
    ..


fi

make clean



