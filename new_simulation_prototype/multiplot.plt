set terminal png size 1920, 1080
set output "_NUM_.png"

set pm3d map
unset xtics
unset ytics
unset key
unset colorbox

#set lmargin at screen 0
#set bmargin  at screen 0
#set tmargin  at screen 1 
#set rmargin  at screen 1

set multiplot


set cbrange[0.8:2.2]
set zrange[0:100]
set xrange[-0.5:0.5]
set yrange[-0.5:0.5]
set size square

set lmargin at screen 0.
set bmargin  at screen 0.5/9.
set tmargin  at screen 8.5/9. 
set rmargin  at screen 0.5 


#_CASE1_ splot "_NUM_.density2D", "_NUM_.faces" u 1:2:(0) w l lc rgb '#FFFFFF' 
#_CASE2_ splot "_NUM_.density2D" 
splot "_NUM_.density2D" 


set cbrange[0.8:2.2]
set zrange[0:100]
set xrange[-0.5:0.5]
set yrange[-0.5:0.5]
set size square

set lmargin at screen 0.5
set bmargin  at screen 0.5/9.
set tmargin  at screen 8.5/9. 
set rmargin  at screen 1. 


#_CASE1_ splot "_NUM_.density2Deuler", "_NUM_.faces" u 1:2:(0) w l lc rgb '#FFFFFF' 
#_CASE2_ splot "_NUM_.density2Deuler" 
splot "_NUM_.density2Deuler" 

