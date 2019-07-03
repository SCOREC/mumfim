#!/bin/bash
  CC=`which mpicc`
  CXX=/gpfs/u/home/PASC/PASCmrsn/scratch/dcs/kokkos/bin/nvcc_wrapper
  #CXX=`which mpicxx`
  cmake \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=0 \
      -DBUILD_TESTS:BOOL=ON \
      -DCMAKE_C_COMPILER=$CC \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_LINKER=`which mpicxx` \
      -DLOGRUN=1 \
      -DENABLE_KOKKOS=ON \
      -DKokkos_DIR=/gpfs/u/home/PASC/PASCmrsn/scratch/dcs/install/kokkos/lib/CMake/Kokkos \
      -DENABLE_WRITE_MICRO_PER_ITER=OFF \
      -DENABLE_WRITE_MICRO_PER_STEP=OFF \
      -DCMAKE_INSTALL_PREFIX=/gpfs/u/home/PASC/PASCmrsn/scratch/dcs/install/biotissue \
      -Dyaml-cpp_DIR=/gpfs/u/home/PASC/PASCmrsn/scratch/dcs/install/yaml-cpp/0.3.0/lib/pkgconfig/ \
      -Damsi_DIR=/gpfs/u/home/PASC/PASCmrsn/scratch/dcs/install/amsi/lib/cmake/amsi/ \
      -DSCOREC_DIR=/gpfs/u/home/PASC/PASCmrsn/scratch/dcs/install/core/lib/cmake/SCOREC/ \
      -Dlas_DIR=/gpfs/u/home/PASC/PASCmrsn/scratch/dcs/install/las/lib/cmake  \
      -Dlas_core_DIR=/gpfs/u/home/PASC/PASCmrsn/scratch/dcs/install/las/lib/cmake \
      -DENABLE_VERBOSITY=OFF \
      -DCMAKE_CXX_STANDARD=11 \
      /gpfs/u/home/PASC/PASCmrsn/barn/biotissue

      #-DCMAKE_C_FLAGS="-O5" \
      #-DCMAKE_CXX_FLAGS="-O5" \
      #-DCMAKE_PREFIX_PATH=/gpfs/u/home/PASC/PASCmrsn/scratch/dcs/install/amsi/lib/cmake/amsi/ \
      #-DCMAKE_CXX_FLAGS="-Wno-unused-variable -Wno-unused-but-set-variable" \
