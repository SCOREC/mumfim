#!/bin/bash

dir=$1
file=$2
prefix=${file%.*}
cd $dir

awk_cmd='{ a=substr($0,1,1); print $0 > ("'
awk_cmd=${awk_cmd}"${prefix}_gantt.\" a \".dat\")}"

awk "$awk_cmd" "$file"

line_count=`verifyLineCount.sh "$dir" "$prefix" "dat"`

file_list=`ls $dir | grep "${prefix}.*.dat"`

action["Active"]="min"
action["Idle"]="max"
o_t=0.0
t=0.0
for (( i = 1; i <= line_count; i++)) ; do
  for f in $file_list; do 
    line=`sed '$iq;d' $f`
    tokens=( $line )
    state=${tokens[3]]}
    if [ "${action[$state]}" == "max" ] ; then
      if (( $(echo "${tokens[2]}" "$t" | awk '{print ($1 > $2)}') )) ; then
        t=${tokens[2]}
      fi
    elif [ "${action[$state]}" == "min" ] ; then
      if (( $(echo "${tokens[2]}" "$t" | awk '{print ($1 < $2)}') )) ; then
        t=${tokens[2]}
      fi
    fi
  done
  echo "0\t${o_t}\t${t}\t${state}"
  o_t=$t
done
