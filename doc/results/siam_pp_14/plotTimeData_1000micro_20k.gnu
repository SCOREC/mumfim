myheight = 400
mywidth = 640
set terminal wxt enhanced size mywidth,myheight
set termoption dash

set title '1st iteration of multiscale solver - 20k mesh - 800 micro processes'

set ylabel 'Time (s)'
set xlabel '# of processes on macroscale'

set grid
set key right box

set logscale xy
set xrange [1.5:38]
set yrange [0.01:8000]

plot \
'allTimes_800micro_20k.dat' every ::1 using ($1):($8) with linespoints title 'Total Solve', \
'allTimes_800micro_20k.dat' every ::1 using ($1):($7) with linespoints title 'Matrix Solve', \
'allTimes_800micro_20k.dat' every ::1 using ($1):($6) with linespoints title 'Assemble', \
'allTimes_800micro_20k.dat' every ::1 using ($1):($3+$5) with linespoints title 'Send+Receive', \
'allTimes_800micro_20k.dat' every ::1 using ($1):($4) with linespoints title 'RVEs'

#'time5188.dat' using (log($1)):(log($2)) with linespoints notitle

pause -1

