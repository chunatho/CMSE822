#!/usr/bin/gnuplot

set term postscript eps enhanced color;
set output '|ps2pdf - intelStream.pdf';
set title "STREAM Bandwidth comparison of Intel compiler settings";
set xlabel "N"; set ylabel "Bytes/sec";
set logscale x;
plot 'O1' title "O1" with lines , \
     'O3' title "O3" with lines, \
     'SSE' title "SSE" with lines, \
     'AVX' title "AVX" with lines, \
     'fast' title "fast" with lines

