#!/bin/bash

dir=$1
prefix=$2
suffix=$3

o_line_count=-1
line_count=-1
for f in `ls $dir | grep "${prefix}.*.${suffix}"` ; do
  line_count=`wc -l < $f`
  if [[ "$o_line_count" -eq "-1" ]] ; then
    o_line_count=line_count
  elif [[ "$o_line_count" -ne "$line_count" ]] ; then
    echo "Error: file lengths are not identical!"
  fi
done

echo $line_count
