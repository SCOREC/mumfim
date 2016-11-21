#!/bin/bash

if [ $# -ne 2 ] ; then
  echo "Usage : plot_weights.sh [directory] [log_file_prefix]"
fi

dir=$1
prefix=$2

cd $dir

# put all the weightings in the same file
cat ${prefix}.*.log > allweights.txt

# split into seperate files based on the first character on each line - the load step number
awk_cmd='{ a=substr($1,0,length($1)-1); print $0 > ("'
awk_cmd=${awk_cmd}"${prefix}_step.\" a \".dat\")}"

awk "$awk_cmd" allweights.txt
rm allweights.txt

#create a plot script and run it
echo "set terminal postscript eps color solid
set output '${prefix}.eps'
unset xtics
set bars fullwidth
set style fill transparent solid 0.25 border -1
set yrange [0:*]
set multiplot layout 2,5 rowsfirst" > $prefix.gnu
for file in `ls ${prefix}_step.*.dat` ; do 
  echo "plot '$file' u 2 w boxes lc rgb\"green\" notitle" >> $prefix.gnu
done
echo "unset multiplot" >> $prefix.gnu



#line_count=`verifyLineCount.sh "${dir}" "${prefix}" "log"`
#for (( i = 2 ; i <= line_count; ii++ )) ; do
#  echo 
#  for file in `ls $dir | grep "${prefix}.*.log"` ; do
#    IFS='.' read -a filename_tokens <<< "$file"
#    process=${filename_tokens[1]}
#    line=`sed '$iq;d' $file`
#    tokens=( $line )
#    echo "${process}\t${tokens[1]}" >> "${prefix}_step_$(( i - 2 )).dat"
#  done
#done
