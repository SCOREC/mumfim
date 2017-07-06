#!/bin/bash -x

# load modules to get packages that biotissue/amsi were compiled with
source $DEVROOT/config_modules.sh

#cdash output root
cd $DEVROOT/test/biotissue/cdash
#remove compilation directories created by previous nightly.cmake runs
rm -rf build/

#run nightly.cmake script
ctest --output-on-failure --script $DEVROOT/test/biotissue/src/core/cdash/jenga-nightly.cmake &> cmake_log.txt

# we will add this in when we start publishing the documentation to the web
# cp cmake_log.txt /net/web/public/seol/scorec/cdash/nightly_cmake_log.txt

#if [ -d "/fasttmp/mersoj/biotissue-test/cdash/build/master" ]; then
#  #core repository checked out by nightly.cmake
#  cd /fasttmp/mersoj/biotissue-test/build/master
#  #build the Doxygen html documentation
#  make doc
#  if [ -d "$PWD/doc/html" ]; then
#    #remove the old web documentation
#    rm -rf /net/web/public/seol/scorec/doxygen
#    #replace it with the generated one
#    cp -r doc/html /net/web/public/seol/scorec/doxygen
#  fi
#fi

# we will add this in when we start using coverty analysis
##core repository checked out by nightly.cmake
#cd /fasttmp/mersoj/biotissue-test/cdash/build/master
##clean the build of object files
#make clean
##run Coverity static analysis on the build
#export PATH=$PATH:/fasttmp/seol/scorec/cov-analysis-linux64-7.7.0.4/bin
#cov-build --dir cov-int make -j 4
##pack up the tarball of results
#tar czvf pumi.tgz cov-int
##cleanup the Chef test output
#cd /fasttmp/seol/scorec
#find meshes/phasta -name "*procs_case" | xargs rm -rf
#find meshes/phasta -name "out_mesh" | xargs rm -rf
