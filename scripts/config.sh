#!/bin/bash
# Cmake config for bio
# usage ./config.sh [build_type] [logrun_flag]
source ./util
ROOT=`dirname $PWD`
LOGRUN_OVERRIDE=$2
if [ -z $1 ]; then
  BUILD_TYPE=Debug
else
  BUILD_TYPE=$1
fi
if [ "$BUILD_TYPE" == "Debug" ] ; then
  BUILD_DIR=$ROOT/build_debug
  BUILD_TESTS=ON
  VERBOSITY_LEVEL=HIGH
elif [ "$BUILD_TYPE" == "Release" ] ; then
  #BUILD_TYPE="RelWithDebugInfo"
  BUILD_DIR=$ROOT/build_release
  BUILD_TESTS=OFF
  VERBOSITY_LEVEL=OFF
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
    -DCMAKE_INSTALL_PREFIX=/gpfs/u/scratch/PASC/shared/install/bio/ \
    -DBUILD_TESTS=OFF \
    -DLOGRUN=$LOGRUN \
    -DCMAKE_PREFIX_PATH=/gpfs/u/scratch/PASC/shared/install/amsi/lib/cmake/amsi \
    -DSCOREC_DIR=/gpfs/u/scratch/PASC/shared/install/core/lib/cmake/SCOREC \
    -DSPARSKIT_DIR=/gpfs/u/scratch/PASC/shared/install/sparskit \
    -DVERBOSITY=0 \
    ..
  chmod g+rw $BUILD_DIR
else
  CC=`which mpicc`
  CXX=/lore/mersoj/kokkos/kokkos/bin/nvcc_wrapper
  #CXX=`which mpicxx`
  cmake \
      -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
      -DBUILD_TESTS:BOOL=ON \
      -DCMAKE_C_COMPILER=$CC \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DLOGRUN=1 \
      -DENABLE_WRITE_MICRO_PER_ITER=OFF \
      -DENABLE_WRITE_MICRO_PER_STEP=OFF \
      -DCMAKE_INSTALL_PREFIX=$DEVROOT/install/bio/$BUILD_TYPE \
      -DCMAKE_PREFIX_PATH="$DEVROOT/install/amsi/$BUILD_TYPE/lib/cmake/amsi/" \
      -DSCOREC_DIR=$DEVROOT/install/core/$BUILD_TYPE/lib/cmake/SCOREC \
      -DKokkos_DIR="/lore/mersoj/kokkos/install/lib/CMake/Kokkos" \
      -Dlas_DIR=$DEVROOT/install/las/RelWithDebugInfo/lib/cmake \
      -Dlas_core_DIR=$DEVROOT/install/las/RelWithDebugInfo/lib/cmake \
      -DENABLE_VERBOSITY=$VERBOSITY_LEVEL \
      -DCMAKE_CXX_STANDARD=11 \
      -DCMAKE_CXX_FLAGS="-Wno-unused-variable -Wno-unused-but-set-variable" \
      ..
fi
      # -DMEMORYCHECK_SUPPRESSIONS_FILE=$DEVROOT/install/openmpi/1.10.7/share/openmpi/openmpi-valgrind.supp \

      #-DCMAKE_CXX_FLAGS="-pg -Wno-unused-variable -Wno-unused-but-set-variable" \
      #-DCMAKE_EXE_LINKER_FLAGS="-pg" \


