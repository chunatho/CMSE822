#!/usr/bin/gnuplot

set term postscript eps enhanced color;
set output '|ps2pdf - intelCompiler_bandwidth.pdf';
set title "Bandwidth comparison of Intel compiler settings";
set xlabel "N"; set ylabel "Bytes/sec";
set logscale x;
plot 'vectO1_output' using 1:4 title "O1" with lines , \
     'vectO3_output' using 1:4 title "O3" with lines, \
     'vectO3SSE_output' using 1:4 title "SSE" with lines, \
     'vectO3AVX_output' using 1:4 title "AVX" with lines, \
     'vectO3fast_output' using 1:4 title "fast" with lines

