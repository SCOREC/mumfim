#!/bin/bash -x

# ./partition $DEVROOT/biotissue/macro/test/dogBone/models/mixed/dogBone_n_bc.smd $DEVROOT/biotissue/macro/test/dogBone/meshes/ 2k db 2

module unload openmpi
module load openmpi/1.6.5-ib
module load simmetrix/simModSuite

MODEL=$1
MESH_ROOT=$2
MESH_SIZE=$3
MESH_NAME=$4
NUM_PARTS=$5

NUM_PROCS=$((NUM_PARTS / 128))
if [ $NUM_PROCS == 0 ]; then
  NUM_PROCS=1
fi

MESH_NAME=${MESH_NAME}_${MESH_SIZE}
if [ $NUM_PROCS == 1 ]; then
  MESH_NAME=${MESH_NAME}.sms
else
  MESH_NAME=${MESH_NAME}_${NUM_PROCS}_parts.sms
fi

mpirun -np $NUM_PROCS \
    bin/x64_rhel5_gcc41/partitionMesh \
    $MODEL \
    $MESH_ROOT/$MESH_SIZE/ \
    $MESH_NAME \
    $NUM_PARTS

