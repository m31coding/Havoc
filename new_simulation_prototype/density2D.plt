set terminal X11
set pm3d map
set size ratio 1
set cbrange[0:6]
set zrange[0:100]
set xrange[-0.5:0.5]
set yrange[-0.5:0.5]
set size square
unset key 
_CASE1_ splot "_NUM_.density2D", "_NUM_.faces" u 1:2:(0) w l lc rgb '#000000' 
_CASE3_ splot "_NUM_.density2D", "_NUM_.faces" u 1:2:(0) w l lc rgb '#000000', "_NUM_.ghost_faces" u 1:2:(0) w l lc rgb '#000000'
_CASE2_ splot "_NUM_.density2D" 
pause -1 "hit return to continue"
