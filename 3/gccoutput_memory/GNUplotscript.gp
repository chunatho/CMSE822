#!/usr/bin/gnuplot

set term postscript eps enhanced color;
set output '|ps2pdf - gccCompiler_bandwidth.pdf';
set title "Bandwidth comparison of GCC compiler settings";
set xlabel "N"; set ylabel "Bytes/sec";
set logscale x;
plot 'vectO1_output' using 1:4 title "O1" with lines, \
     'vectO3_output' using 1:4 title "O3" with lines, \
     'vectO3march_output' using 1:4 title "O3 march" with lines

