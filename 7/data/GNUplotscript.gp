#!/usr/bin/gnuplot

set term postscript eps enhanced color 
set title "Time to complete 100 temporal iterations";
set output '|ps2pdf - versioncomp.pdf';
set logscale y; set logscale x;
set xlabel "length of 2D array (n)"; set ylabel "Time (sec)";

plot 'v1times' using 1:2 title "seperate kernels" lt rgb "red", \
     'v2times' using 1:2 title "colalesced kernels" lt rgb "black", \
     'v3times' using 1:2 title "tmp variables" lt rgb "blue"

set term postscript eps enhanced color 
set title "Comparison of accelarted code";
set output '|ps2pdf - accelcomp.pdf';
set logscale y; set logscale x;
set xlabel "length of 2D array (n)"; set ylabel "Time (sec)";

plot 'v0times' using 1:2 title "non-accel" lt rgb "red", \
     'v3times' using 1:2 title "CUDA" lt rgb "black"


