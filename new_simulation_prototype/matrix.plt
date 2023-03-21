set terminal postscript eps color size 6.2,6.2
set output "_NUM_.eps"

#set pm3d map interpolate 0,0
#set pm3d map
set size ratio 1
set zrange[0:100]
set cbrange[0:6]
unset xtics
unset ytics
unset xlabel
unset ylabel
unset title 
unset key

plot "_NUM_.density_sample" matrix with image
#pause -1 "hit return"
