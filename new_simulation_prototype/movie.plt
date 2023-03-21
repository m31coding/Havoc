set terminal png size 1280,1024
set output "_NUMBER_.png"

h0 = 0.05
w0 = 0.12
h = 0.9
w = 1024.0/1280.0*h

set lmargin at screen w0
set rmargin at screen w0 + w

set bmargin at screen h0
set tmargin at screen h0 + h


set title "time: _TIME_s (_NUMBER_)"
set pm3d map
set size ratio 1
set xrange[-0.5:0.5]
set yrange[-0.5:0.5]
unset key
set cblabel "surface density"
set cbrange[0:6]
set zrange[0:100]
splot "_GRAPH_", "_FACES_" u 1:2:(0) w l lc rgb '#FFFFFF' 
#splot "_GRAPH_" 

