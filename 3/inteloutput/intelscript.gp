#!/usr/bin/gnuplot

set term postscript eps enhanced color;
set output '|ps2pdf - intelCompiler.pdf';
set title "Performance comparison of Intel compiler settings";
set xlabel "N"; set ylabel "GFLOPS/sec";
set logscale x;
plot 'vectO1_output' using 1:4 title "O1", \
     'vectO3_output' using 1:4 title "O3", \
     'vectO3SSE_output' using 1:4 title "SSE", \
     'vectO3AVX_output' using 1:4 title "AVX", \
     'vectO3fast_output' using 1:4 title "fast"

