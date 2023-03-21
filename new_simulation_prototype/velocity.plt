set xrange[0:1]
set yrange[-0.05:0.05]

scale = _SCALE_ 
plot "_FILE_" u 1:2:($3*scale):($4*scale) w vectors
pause -1
