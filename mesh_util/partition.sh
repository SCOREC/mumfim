#!/bin/bash -x

if [[ $# < 5 ]] ; then
  echo "Usage : $0 [.smd] [dir] [mesh_suffix] [prefix] [num_parts] {debug}"
  exit
fi
# ./partition.sh $DEVROOT/biotissue/macro/test/dogBone/models/mixed/dogBone_n_bc.smd $DEVROOT/biotissue/macro/test/dogBone/meshes 2k db 2

MODEL=$1
MESH_ROOT=$2
MESH_SIZE=$3
MESH_NAME=$4
NUM_PARTS=$5
DEBUG=$6

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

if [ "$DEBUG" == "Debug" ] ; then
  TV="-tv"
fi

mpirun -np $NUM_PROCS $TV \
    bin/x64_rhel7_gcc48/partitionMesh \
    $MODEL \
    $MESH_ROOT/$MESH_SIZE \
    $MESH_NAME \
    $NUM_PARTS

