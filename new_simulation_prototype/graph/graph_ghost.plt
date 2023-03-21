set terminal X11
grey = "lc rgb '#808080'"
green = "lc rgb '#228B22'"
blue = "lc rgb '#4169E1'"
orchid = "lc rgb '#DA70D6'"
orangered = "lc rgb '#FF4500'"
red = "lc rgb '#FF0000'"
darkblue = "lc rgb '#191970'"
black = "lc rgb '#000000'"
brown = "lc rgb '#A0522D'"
eps = "set terminal postscript eps color enhanced font ',20'"
map = "set pm3d map; set palette rgb 33,13,10; set border 10"
set size ratio 1 
plot "ghost_faces" w l lc rgb '#C9C9C9', "faces" w l lc rgb '#228B22', "sites" w p lc rgb '#FF0000' pt 2, "ghosts" w p lc rgb '#F9966B' pt 2, "centers" w p @black pt 4, "ghost_centers" w p pt 4 lc rgb '#C9C9C9', "vertices" w p lc rgb '#4169E1' pt 2, "breakpoint_vertices" w p lc rgb '#9778E1' pt 2, "breakpoint_edges" w l lc rgb '#00c44B', "face_centers" w p lc rgb '#228B22' pt 6
pause -1 "hit return to continue"
