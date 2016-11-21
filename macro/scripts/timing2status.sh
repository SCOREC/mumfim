#!/bin/bash                                                                                      

#1st argument is input file name (without .#.log)
#2nd argument is number of microprocessors
#3rd argument is output file name (without .#.log)

for ((i=0;i<$2;i++)) 
do
    ./timing2gantt.sh $1.$i.log > $3.$i.log
done
