#!/usr/bin/gnuplot

set term postscript eps enhanced color;
set output '|ps2pdf - gccStream.pdf';
set title "STREAM Bandwidth comparison of GCC compiler settings";
set xlabel "N"; set ylabel "Bytes/sec";
set logscale x;
plot 'gccO1' title "O1" with lines, \
     'gccO3' title "O3" with lines, \
     'gccmarch' title "march" with lines

