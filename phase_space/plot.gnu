#set term postscript eps enhanced color "Helvetica" 36
set term postscript eps solid enhanced color "Helvetica" 36
#set term postscript eps solid enhanced color "Times" 36
set output "plot.eps"

set style line 1 lt rgb "red" lw 7 # red
#set style line 2 lt rgb "red" lw 7 dt 2 # red dashed
set style line 3 lt rgb "green" lw 7 # green
#set style line 4 lt rgb "green" lw 7 dt 2 # green dashed
set style line 5 lt rgb "blue" lw 7 # blue
#set style line 6 lt rgb "blue" lw 7 dt 2 # blue dashed
set style line 7 lt rgb "black" lw 7 # black
#set style line 8 lt rgb "black" lw 7 dt 2 # black dashed
set style line 9 lt rgb "grey" lw 7 # grey
#set style line 9 lt rgb "grey" lw 7 dt 2 # grey dashed
set style line 10 lt rgb "black" lw 1 # black
set style line 11 lt rgb "red" lw 1 # black

set key top right
#set label at 0.5,5.75 "{/Symbol a}=4"

set multiplot

set size ratio 0.6
#set lmargin 3
#set bmargin 3
#set rmargin 3
#set tmargin 0
set pointsize 2

set origin 0.0,0.0
set size square
set grid
#set xrange[0:58]
set yrange[-200:200]
#set ytics (5e16,1e17,2e17,3e17,4e17)
set ytics 100
set xtics 2
#unset xtics
#set format y "%1.0t.10^{%L}"
#set format y "10^{%L}"
set ylabel "y'(mrad)" offset 0,0.
#set log x

set xlabel "y(mm)" offset 0,0.25

pl "< awk '$4==3' ../DATA/beam_phase_space.dat" u ($1*1e3):(atan($3/$2)*1e3) w d ls 11 noti, "< awk '$4==4' ../DATA/beam_phase_space.dat" u ($1*1e3):(atan($3/$2)*1e3) w d ls 10 noti 
