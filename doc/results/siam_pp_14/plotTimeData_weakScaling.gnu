myheight = 800
mywidth = 1280
set terminal wxt enhanced size mywidth,myheight linewidth 5 font ",24"
set termoption dash

set title '1st iteration of multiscale solver - Weak Scaling'

reset
set xtics nomirror
set x2tics
set autoscale xfix
set autoscale x2fix

set ylabel 'Time (s)'
set xlabel '# of processes on microscale'
set x2label '# of processes on macroscale'

set grid
set key bottom right box

set logscale xy
set logscale x2
#set xrange [100:10000]
#set yrange [0.008:10000]

set out 'weakScaling.png'

plot \
'allTimes_weakScaling.dat' every ::1 using ($2):($8) with linespoints title 'Total Solve', \
'allTimes_weakScaling.dat' every ::1 using ($2):($7) with linespoints title 'Matrix Solve', \
'allTimes_weakScaling.dat' every ::1 using ($2):($6) with linespoints title 'Assemble', \
'allTimes_weakScaling.dat' every ::1 using ($2):($3+$5) with linespoints title 'Send+Receive', \
'allTimes_weakScaling.dat' every ::1 using ($2):($4) with linespoints title 'RVEs', \
  '' using ($2):($4):x2tic(1) axes x2y1 lt 5 pt 5 notitle

#'time5188.dat' using (log($1)):(log($2)) with linespoints notitle

pause -1

