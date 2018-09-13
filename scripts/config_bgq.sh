#!/bin/bash
  cmake \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=0 \
      -DBUILD_TESTS:BOOL=OFF \
      -DCMAKE_C_COMPILER=`which mpicc` \
      -DCMAKE_CXX_COMPILER=`which mpicxx` \
      -DLOGRUN=1 \
      -DENABLE_WRITE_MICRO_PER_ITER=OFF \
      -DENABLE_WRITE_MICRO_PER_STEP=OFF \
      -DCMAKE_INSTALL_PREFIX=/gpfs/u/home/PASC/PASCmrsn/scratch/install/biotissue \
      -DCMAKE_PREFIX_PATH="/gpfs/u/home/PASC/PASCmrsn/scratch/install/amsi/lib/cmake/amsi/" \
      -DSCOREC_DIR=/gpfs/u/home/PASC/PASCmrsn/scratch/install/core/lib/cmake/SCOREC/ \
      -Dlas_DIR=/gpfs/u/home/PASC/PASCmrsn/scratch/install/las/lib/cmake \
      -Dlas_core_DIR=/gpfs/u/home/PASC/PASCmrsn/scratch/install/las/lib/cmake \
      -DENABLE_VERBOSITY=MED \
      /gpfs/u/home/PASC/PASCmrsn/barn/biotissue
