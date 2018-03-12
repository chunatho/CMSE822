#!/usr/bin/gnuplot

set term postscript eps enhanced color 
set title "Subset performance with respect to Number of proccesses";
set output '|ps2pdf - subset.pdf';
set logscale y;
set xlabel "number of processes"; set ylabel "completetion time (sec)";

plot 'subsetdata_20' title "20 Elements" lt rgb "red", \
     'subsetdata_24' title "24 Elements" lt rgb "black", \
     'subsetdata_28' title "28 Elements" lt rgb "blue"



set term postscript eps enhanced color 
set title "Comparison of block and cyclic decoders";
set output '|ps2pdf - decoder.pdf';
set logscale y;
set xlabel "number of processes"; set ylabel "completetion time (sec)";

plot 'codedata_cycle' title "cycle" lt rgb "red", \
     'codedata_block' title "block" lt rgb "black", \
     'codedata_bonus' title "OpenMP" lt rgb "blue"
