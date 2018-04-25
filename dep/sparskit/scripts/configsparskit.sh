#!/bin/bash
# CMake build for sparskit_scorec
#
# usage: configsparskit.sh [build_type] [install_prefix]
#
# todo: verify that $2 is a valid install prefix
#

#module load las_0.1

HOSTNAME=`hostname`

if [ "$1" == "Debug" -o "$1" == "Release" ]; then

  cd ../build

  cmake \
    -DCMAKE_C_COMPILER=mpicc\
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DCMAKE_Fortran_COMPILER=mpif90\
    -DCMAKE_BUILD_TYPE=$1 \
    -DCMAKE_INSTALL_PREFIX=$2 \
    -DLAS_INCLUDE_DIR=$LAS_INCLUDE_DIR \
    ..


fi

make clean




