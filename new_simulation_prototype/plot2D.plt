set pm3d map
set size ratio 1
set cbrange[0:2]
set zrange[0:2]
splot "_GRAPH_", "_FACES_" u 1:2:(0) w l lc rgb '#FFFFFF' 
pause -1 "hit return to continue"
