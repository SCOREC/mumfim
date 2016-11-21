myheight = 400
mywidth = 640
set terminal wxt enhanced size mywidth,myheight
set termoption dash

set title '1st iteration of multiscale solver - 200k mesh - 16 macro processes'

set ylabel 'Time (s)'
set xlabel '# of processes on microscale'

set grid
set key bottom right box

set logscale xy
set xrange [450:18000]
set yrange [0.02:10000]

plot \
'allTimes_200k.dat' every ::8::13 using ($2):($8) with linespoints title 'Total Solve', \
'allTimes_200k.dat' every ::8::13 using ($2):($7) with linespoints title 'Matrix Solve', \
'allTimes_200k.dat' every ::8::13 using ($2):($6) with linespoints title 'Assemble', \
'allTimes_200k.dat' every ::8::13 using ($2):($3+$5) with linespoints title 'Send+Receive', \
'allTimes_200k.dat' every ::8::13 using ($2):($4) with linespoints title 'RVEs'

#'time5188.dat' using (log($1)):(log($2)) with linespoints notitle

pause -1

