h0 = 0.05
w0 = 0.12
h = 0.9
w = 1024.0/1280.0*h

set lmargin at screen w0
set rmargin at screen w0 + w

set bmargin at screen h0
set tmargin at screen h0 + h

set size ratio 1
set xrange[-0.5:0.5]
set yrange[-0.5:0.5]
set xtics 0.1
set ytics 0.1
set size square
unset grid
unset key
set title "Voronoi particles - _NUM_"

_CASE1_set terminal png size 1280,1024
_CASE1_set output "_NUMBER_particles.png"

plot "_NUM_.data" u 1:2 w d 
#plot "_NUM_.data" u 1:2 w d, "_NUM_.ghosts" u 1:2 w d @blue

_CASE2_pause -1 "hit return to continue"

