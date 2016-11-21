myheight = 400
mywidth = 640
set terminal wxt enhanced size mywidth,myheight
set termoption dash

set title '1st iteration of multiscale solver - 200k mesh - Total solve time'

set ylabel 'Time (s)'
set xlabel '# of processes on microscale'

set grid
set key bottom right box

set logscale xy
set xrange [450:18000]
#set yrange [0.008:10000]

plot \
'allTimes_200k.dat' every ::1::5 using ($2):($8) with linespoints title '4 macro', \
'allTimes_200k.dat' every ::6::9 using ($2):($8) with linespoints title '8 macro', \
'allTimes_200k.dat' every ::10::15 using ($2):($8) with linespoints title '16 macro', \
'allTimes_200k.dat' every ::16::20 using ($2):($8) with linespoints title '32 macro', \
'allTimes_200k.dat' every ::21::23 using ($2):($8) with linespoints title '64 macro'

#'time5188.dat' using (log($1)):(log($2)) with linespoints notitle

pause -1

