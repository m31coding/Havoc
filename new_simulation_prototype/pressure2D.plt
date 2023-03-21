set terminal X11
set pm3d map
set size ratio 1
set cbrange[0:5]
set zrange[0:5]
set xrange[-0.5:0.5]
set yrange[-0.5:0.5]
set size square
_CASE1_ splot "_GRAPH_", "_FACES_" u 1:2:(0) w l lc rgb '#000000' 
_CASE2_ splot "_GRAPH_" 
pause -1 "hit return to continue"
