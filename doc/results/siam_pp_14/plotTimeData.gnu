myheight = 400
mywidth = 640
set terminal wxt enhanced size mywidth,myheight
set termoption dash

set title '1st iteration of multiscale solver - 5k mesh'

set ylabel 'Time (s)'
set xlabel '# of processes on microscale'

set grid
set key right box

set logscale xy
set xrange [10:4300]

plot \
'allTimes_1macro_5k.dat' every ::1 using ($2):($6) with linespoints title 'Total Solve', \
'allTimes_1macro_5k.dat' every ::1 using ($2):($5) with linespoints title 'Matrix Solve', \
'allTimes_1macro_5k.dat' every ::1 using ($2):($4) with linespoints title 'Assemble', \
'allTimes_1macro_5k.dat' every ::1 using ($2):($3) with linespoints title 'Send'

#'time5188.dat' using (log($1)):(log($2)) with linespoints notitle

pause -1

