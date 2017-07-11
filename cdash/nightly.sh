#!/bin/bash

# load modules to get packages that biotissue/amsi were compiled with
source $DEVROOT/config_modules.sh

#cdash output root
cd $DEVROOT/test/biotissue/cdash
#remove compilation directories created by previous nightly.cmake runs
rm -rf build/

#run nightly.cmake script
ctest --output-on-failure --script $DEVROOT/biotissue/cdash/jenga-nightly.cmake &> cmake_log.txt
