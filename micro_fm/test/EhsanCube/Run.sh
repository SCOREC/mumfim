#!/bin/sh -x
$MPI_HOME/bin/mpirun -np 2 -tv\
         $DEVROOT/micro_fm/test/SourceMain/MICRO_FM \
	-s $DEVROOT/micro_fm/test/EhsanCube/cube_bc.smd \
	-m $DEVROOT/micro_fm/test/EhsanCube/mesh_2_parts.sms \
	-l /net/common/meshSim/license/license.txt
