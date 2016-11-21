#!/bin/bash

if [ $# -ne 2 ] ; then
  echo "Usage : ./timing2fig.sh [directory] [log_file_prefix]"
  exit
fi

dir=$1
name=$2

full=$dir'/'$name

$PWD/generate_gantt.sh $dir $name > $full'.dat'
python $PWD/gantt.py -o $full'.gpl' $full'.dat'

echo 'set terminal postscript eps color solid' > $full.gnu
echo 'set output "'${full}.eps'"' >> $full.gnu
echo 'load "'${full}.gpl'"' >> $full.gnu
echo 'unset output' >> $full.gnu

gnuplot $full.gnu
