#!/bin/bash                                                                                      

dir=/fasttmp/wtobin/develop/results/dogbone_20k_11
macro=4
micro="64 128 256 512"

for i in $micro 
do
    echo ${dir}/${macro}_${i}/micro_fo_efficiency
    echo ${dir}/${macro}_${i}/macro_efficiency

    python plot_gantt_v3.py ${dir}/${macro}_${i}/macro_efficiency ${macro} ${dir}/${macro}_${i}/micro_fo_efficiency ${i} average gantt_${macro}\_${i}\_avg_no_migration

    python plot_gantt_v3.py ${dir}/${macro}_${i}/macro_efficiency ${macro} ${dir}/${macro}_${i}/micro_fo_efficiency ${i} max_min gantt_${macro}\_${i}\_maxmin_no_migration   
done
