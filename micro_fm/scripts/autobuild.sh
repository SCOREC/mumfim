#!/bin/bash
cd ../build

module load petsc
module load apf_1.0

dep_libs="$PETSC_DIR/$PETSC_ARCH/lib/libsuperlu_4.3.a"
dep_libs="$dep_libs $PETSC_DIR/$PETSC_ARCH/lib/libsuperlu_dist_3.3.a"
dep_libs="$dep_libs $PETSC_DIR/$PETSC_ARCH/lib/libmetis.so"
dep_libs="$dep_libs $PETSC_DIR/$PETSC_ARCH/lib/libparmetis.so"
dep_libs="$dep_libs $PETSC_DIR/$PETSC_ARCH/lib/libpetsc.so"
dep_libs="$dep_libs $APF_LIB_DIR/libapf.a"
dep_libs="$dep_libs $APF_LIB_DIR/libapf_sim.a"
dep_libs="$dep_libs $DEVROOT/install/lib/libAMSI.a"
dep_libs="$dep_libs $DEVROOT/install/lib/libSCORECUtil.a"
dep_libs="$dep_libs $DEVROOT/install/lib/libFEA.a"

notifybuild \
../src \
"${dep_libs}" \
'make install' \
MICRO_FM
