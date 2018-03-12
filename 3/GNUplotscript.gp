#!/usr/bin/gnuplot

set term postscript eps enhanced color;
set output '|ps2pdf - gccCompiler.pdf';
set title "Performance comparison of GCC compiler settings"; set xlabel "N"; set ylabel "GFLOPS/sec";
set logscale x;
plot 'vectO1_output' using 1:4 title "O1", \
     'vectO3_output' using 1:4 title "O3", \
     'vectO3march_output' using 1:4 title "O3 march"

