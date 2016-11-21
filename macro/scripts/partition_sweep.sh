#!/bin/bash -x

# [model] [mesh_directory] [input_mesh] [num_parts]

for SIZE in 2 4 8 16 32 64 128
do
  partitionMesh $1 $2 $3 $SIZE
done

MESH=${3%????}

for SIZE in 256 512
do
  DIV=$(( $SIZE / 128 ))
  mpirun -np $DIV partitionMesh $1 $2 ${MESH}_${DIV}_parts.sms $SIZE
done
