set terminal png size 1920,1080
set output "_NUMBER_.png"

set lmargin at screen 0
set bmargin  at screen 0
set tmargin  at screen 1 
set rmargin  at screen 1

#set title "time: _TIME_s (_NUMBER_)"
unset border
set pm3d map
set size ratio 1
set xrange[-0.5:0.5]
set yrange[-0.5:0.5]
unset xtics
unset ytics
unset key
#set cblabel "surface density"
unset colorbox
set cbrange[0:6]
set zrange[0:100]
splot "_GRAPH_", "_FACES_" u 1:2:(0) w l lw 1.8 lc rgb '#FFFFFF' 
#splot "_GRAPH_" 

