#!/bin/bash
# CMake build for Sparskit as a Biotissue_AMSI dependency

if [ ! -d ../build ]; then
  mkdir ../build
fi

cd ../build

cmake \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_INSTALL_PREFIX=$DEVROOT/install \
  ..



