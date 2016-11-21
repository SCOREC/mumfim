#!/bin/bash
if [ -z $1 ]; then
  BUILD_TYPE=Debug
else
  BUILD_TYPE=$1
fi

AMSI_LIB=${DEVROOT}/install/amsi/openmpi-1.10.0/lib/libcontrol

if [ "$BUILD_TYPE" == "Debug" ]; then
  BUILD_DIR=../build_debug
  AMSI_LIB=${AMSI_LIB}d
elif [ "$BUILD_TYPE" == "Release" ] ; then
  BUILD_DIR=../build_release
fi

AMSI_LIB=${AMSI_LIB}.a

if [ ! -d $BUILD_DIR ]; then
  mkdir $BUILD_DIR
fi
cd $BUILD_DIR

notifybuild \
"../macro/src/ ../micro_fo/src/ ../micro_fm/src/ ../macro/test/multiscale/ ../macro/test/singlescale ../micro_fo/test/" \
"${AMSI_LIB}" \
'make' \
"BIOTISSUE_$BUILD_TYPE"

