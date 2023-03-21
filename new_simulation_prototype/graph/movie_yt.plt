set terminal png size 1080,1080
set output "_NUMBER_.png"

set lmargin at screen 0
set bmargin  at screen 0
set tmargin  at screen 1 
set rmargin  at screen 1

line_width=2
point_size=0.8
set size ratio 1

#point line width
pointlw = 2.5

unset grid
unset xtics
unset ytics
unset key
unset border
set size square

set xrange[-0.68:-0.3]
set yrange[0.3:0.68]

#set object 1 rectangle from -0.62,-0.62 to 0.62, 0.62 behind lw 5
#set object 2 rectangle from -0.5,-0.5 to 0.5,0.5 behind lt 0 

#real box
set arrow from -0.5,0.3 to -0.5,0.5 nohead lt 1 lw 5 lc rgb "green" 
set arrow from -0.5,0.5 to -0.3,0.5 nohead lt 1 lw 5 lc rgb "green" 

#ghost box
set arrow from -0.62,0.3 to -0.62,0.62 nohead lt 1 lw 5 lc rgb "green" 
set arrow from -0.62,0.62 to -0.3,0.62 nohead lt 1 lw 5 lc rgb "green" 

plot "ghosts" w p lc @color10 ps point_size lw pointlw pt 2,\
"ghosts_fixed" w p lc @color14 ps point_size lw pointlw pt 2,\
"sites" w p lc @color4 ps point_size lw pointlw pt 2,\
"faces" w l lc rgb '#000000' lw line_width lt 1,\
"ghost_faces" w l lc rgb '#000000' lw line_width lt 1,\
"breakpoint_edges" w l lc rgb '#00000' lw line_width lt 1
