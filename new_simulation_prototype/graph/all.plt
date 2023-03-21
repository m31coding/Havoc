unset clip points
set clip one
unset clip two
set bar 1.000000 front
set border 3 back linetype 1 linecolor rgb "#808080"  linewidth 1.500
set timefmt z "%d/%m/%y,%H:%M"
set zdata 
set timefmt y "%d/%m/%y,%H:%M"
set ydata 
set timefmt x "%d/%m/%y,%H:%M"
set xdata 
set timefmt cb "%d/%m/%y,%H:%M"
set timefmt y2 "%d/%m/%y,%H:%M"
set y2data 
set timefmt x2 "%d/%m/%y,%H:%M"
set x2data 
set boxwidth
set style fill  empty border
set style rectangle back fc  lt -3 fillstyle   solid 1.00 border lt -1
set style circle radius graph 0.02, first 0, 0 
set style ellipse size graph 0.05, 0.03, first 0 angle 0 units xy
set dummy x,y
set format x "% g"
set format y "% g"
set format x2 "% g"
set format y2 "% g"
set format z "% g"
set format cb "% g"
set format r "% g"
set angles radians
set tics back
set grid nopolar
set grid xtics nomxtics ytics nomytics noztics nomztics \
 nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics
set grid back   linetype 0 linecolor rgb "#808080"  linewidth 1.000,  linetype 0 linecolor rgb "#808080"  linewidth 1.000
set raxis
set key title ""
set key inside right top vertical Right noreverse enhanced autotitles nobox
set key noinvert samplen 4 spacing 1 width 0 height 0 
set key maxcolumns 0 maxrows 0
set key noopaque
unset key
unset label
unset arrow
set style increment default
unset style line
set style line 1  linetype 1 linewidth 1.500 pointtype 1 pointsize default pointinterval 0
set style line 50  linetype 1 linecolor rgb "#808080"  linewidth 1.500 pointtype 1 pointsize default pointinterval 0
set style line 51  linetype 0 linecolor rgb "#808080"  linewidth 1.000 pointtype 0 pointsize default pointinterval 0
unset style arrow
set style histogram clustered gap 2 title  offset character 0, 0, 0
unset logscale
set offsets 0, 0, 0, 0
set pointsize 1
set pointintervalbox 1
set encoding default
unset polar
unset parametric
unset decimalsign
set view 60, 30, 1, 1
set samples 100, 100
set isosamples 10, 10
set surface
unset contour
set clabel '%8.3g'
set macros
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels auto 5
set cntrparam points 5
set size ratio 1 1,1
set origin 0,0
set style data points
set style function lines
set xzeroaxis linetype -2 linewidth 1.000
set yzeroaxis linetype -2 linewidth 1.000
set zzeroaxis linetype -2 linewidth 1.000
set x2zeroaxis linetype -2 linewidth 1.000
set y2zeroaxis linetype -2 linewidth 1.000
set ticslevel 0.5
set mxtics default
set mytics default
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set xtics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set xtics autofreq  norangelimit
set ytics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set ytics autofreq  norangelimit
set ztics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set ztics autofreq  norangelimit
set nox2tics
set noy2tics
set cbtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set cbtics autofreq  norangelimit
set rtics axis in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set rtics autofreq  norangelimit
set title "" 
set title  offset character 0, 0, 0 font "" norotate
set timestamp bottom 
set timestamp "" 
set timestamp  offset character 0, 0, 0 font "" norotate
set rrange [ * : * ] noreverse nowriteback
set trange [ * : * ] noreverse nowriteback
set urange [ * : * ] noreverse nowriteback
set vrange [ * : * ] noreverse nowriteback
set xlabel "" 
set xlabel  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set x2label "" 
set x2label  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set xrange [ -3.17858 : -1.09642 ] noreverse nowriteback
set x2range [ -47.9118 : -38.3038 ] noreverse nowriteback
set ylabel "" 
set ylabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set y2label "" 
set y2label  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set yrange [ 4.24086 : 6.32301 ] noreverse nowriteback
set y2range [ 7.51068 : 12.2333 ] noreverse nowriteback
set zlabel "" 
set zlabel  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set zrange [ * : * ] noreverse nowriteback
set cblabel "" 
set cblabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set cbrange [ * : * ] noreverse nowriteback
set zero 1e-08
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
set locale "C"
set pm3d explicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles nohidden3d corners2color mean
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 front bdefault
set style boxplot candles range  1.50 outliers pt 7 separation 1 labels auto unsorted
set loadpath 
set fontpath 
set psdir
set fit noerrorvariables
GNUTERM = "x11"
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
plot "faces" w lines lc rgb '#228B22' lw 2, "sites" w p ps 2 lw 2lc rgb '#FF0000', "vertices" w p ps 2 pt 7
