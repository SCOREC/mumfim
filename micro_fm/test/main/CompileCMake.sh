
# CMake build for SourceMain

HOSTNAME=`hostname`

if [ "$HOSTNAME" == "q.ccni.rpi.edu" ]; then
  cmake \
    -DHOSTNAME=bgq \
    -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=/gpfs/u/home/PASC/PASCfvgd/scratch-shared/biotissue/bin \
    -DPETSC_PATH=/gpfs/u/home/PASC/PASCfvgd/barn-shared/petsc-3.2-p7 \
    -DPETSC_ARCH=arch-linux2-c-opt \
    -DBIOTISSUE_LIBRARY_PATH=$DEVROOT/install/lib \
    -DBOOST_INCLUDE_PATH=/gpfs/u/home/PASC/PASCfvgd/barn-shared/boost_1_44_0\
    -DAMSI_INCLUDE_PATH=$DEVROOT/install/include \
    -DMESHSIM_INCLUDE_PATH=/gpfs/u/home/PASC/PASCfvgd/barn-shared/simmetrix/9.0-130517 \
    -DMESHSIM_LIB_SUBDIR=ppc64_rhel6_xlc \
    -DSCORECUTIL_INCLUDE_PATH=$DEVROOT/install/include \
    -DFEMANALYSIS_INCLUDE_PATH=$DEVROOT/install/include \
    -DSPARSKIT_PATH=$DEVROOT/install \
    -MICRO_FM_INCLUDE_PATH=$DEVROOT/install/include \
    -MICRO_FO_INCLUDE_PATH=$DEVROOT/install/include \
    -DMANUAL_SET_MPI=1 \
    -DMPI_INCLUDE_PATH=/bgsys/drivers/V1R2M0/ppc64/comm/gcc/include \
    -DMPI_COMPILER=/bgsys/drivers/V1R2M0/ppc64/comm/gcc/bin/mpic++ \
    -DMPI_LINK_FLAGS=-Wl,--export-dynamic\
    .
else
  cmake \
    -DHOSTNAME=scorec \
    -DBIOTISSUE_LIBRARY_PATH=$DEVROOT/install/lib \
    -DPETSC_PATH=$DEVROOT/dep/petsc-3.2-p7 \
    -DPETSC_ARCH=arch-linux2-c-debug \
    -DBOOST_INCLUDE_PATH=$DEVROOT/dep/boost/boost_1_44_0 \
    -DAMSI_INCLUDE_PATH=$DEVROOT/install/include \
    -DMESHSIM_INCLUDE_PATH=/net/common/meshSim/latest/ \
    -DMESHSIM_LIB_SUBDIR=x64_rhel5_gcc41 \
    -DFEA_INCLUDE_PATH=$DEVROOT/install/include \
    -DMANUAL_SET_MPI=0 \
    .
fi

make clean
make VERBOSE=1




