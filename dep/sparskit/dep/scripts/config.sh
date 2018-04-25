#!/bin/bash
# CMake build for Sparskit as a Biotissue_AMSI dependency

if [ ! -d ../build ]; then
  mkdir ../build
fi

cd ../build

cmake \
  -DCMAKE_C_COMPILER=mpicc\
  -DCMAKE_CXX_COMPILER=mpicxx\
  -DCMAKE_Fortran_COMPILER=mpif90\
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$DEVROOT/install/sparskit/ \
  ..



