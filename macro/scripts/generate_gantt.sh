#!/bin/bash

dir=$1
prefix=$2

for file in `ls $dir | grep "${prefix}.*.log"`
do
  IFS='.' read -a tokens <<< "$file"
  process=${tokens[1]}
  if [ "$process" -le 64 ]; then
    $PWD/timing2gantt.sh $dir'/'$file $process
  fi
done
