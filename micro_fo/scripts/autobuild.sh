#!/bin/bash
cd ../build

module load AMSI_0.1
module load apf_1.0
module load Sparskit_1.0
module load SCORECUTIL_1.0
module load FEA_0.1

dep_libs="$dep_libs $AMSI_LIB_DIR/libAMSI.a"
dep_libs="$dep_libs $APF_LIB_DIR/libapf.a"
dep_libs="$dep_libs $APF_LIB_DIR/libapf_sim.a"
dep_libs="$dep_libs $SPARSKIT_LIB_DIR/libSPARSKIT.a"
dep_libs="$dep_libs $SCORECUTIL_LIB_DIR/libSCORECUtil.a"
dep_libs="$dep_libs $FEA_LIB_DIR/libFEA.a"

notifybuild \
../src \
"${dep_libs}" \
'make install' \
MICRO_FO
