#!/bin/bash
  CC=`which mpicc`
  #CXX=/lore/mersoj/kokkos/kokkos/bin/nvcc_wrapper
  CXX=`which mpicxx`
  cmake \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=0 \
      -DBUILD_TESTS:BOOL=ON \
      -DCMAKE_C_COMPILER=$CC \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DLOGRUN=1 \
      -DENABLE_KOKKOS=OFF \
      -DENABLE_WRITE_MICRO_PER_ITER=OFF \
      -DENABLE_WRITE_MICRO_PER_STEP=OFF \
      -DENABLE_KOKKOS=OFF \
      -DCMAKE_INSTALL_PREFIX=/gpfs/u/home/PASC/PASCmrsn/scratch/install-erp/biotissue \
      -DCMAKE_PREFIX_PATH=/gpfs/u/home/PASC/PASCmrsn/scratch/install-erp/amsi/lib/cmake/amsi/ \
      -DSCOREC_DIR=/gpfs/u/home/PASC/PASCmrsn/scratch/install-erp/core/lib/cmake/SCOREC/ \
      -Dlas_DIR=/gpfs/u/home/PASC/PASCmrsn/scratch/install-erp/las/lib/cmake  \
      -Dlas_core_DIR=/gpfs/u/home/PASC/PASCmrsn/scratch/install-erp/las/lib/cmake \
      -DENABLE_VERBOSITY=$VERBOSITY_LEVEL \
      -DCMAKE_CXX_STANDARD=11 \
      -DCMAKE_CXX_FLAGS="-Wno-unused-variable -Wno-unused-but-set-variable" \
      /gpfs/u/home/PASC/PASCmrsn/barn/biotissue
