#!/usr/bin/gnuplot

set term postscript eps enhanced color 
set title "Comparison of Time to Sort";
set output '|ps2pdf - SORTcomp.pdf';
set logscale y;
set xlabel "number of processes"; set ylabel "Time (sec)";

plot 'v1times_1' using 1:2 title "size1 Single RNG" lt rgb "red", \
     'v2times_1' using 1:2 title "size1 Distributed RNG" lt rgb "black", \
     'v1times_2' using 1:2 title "size2 Single RNG" lt rgb "blue", \
     'v2times_2' using 1:2 title "size2 Distributed RNG" lt rgb "cyan", \
     'v1times_3' using 1:2 title "size3 Single RNG" lt rgb "green", \
     'v2times_3' using 1:2 title "size3 Distributed RNG" lt rgb "violet"


set term postscript eps enhanced color 
set title "Comparison of Time to Generate/Distribute RNGs";
set output '|ps2pdf - RNGcomp.pdf';
set logscale y;
set xlabel "number of processes"; set ylabel "Time (sec)";

plot 'v1times_1' using 1:3 title "size1 Single RNG" lt rgb "red", \
     'v2times_1' using 1:3 title "size1 Distributed RNG" lt rgb "black", \
     'v1times_2' using 1:3 title "size2 Single RNG" lt rgb "blue", \
     'v2times_2' using 1:3 title "size2 Distributed RNG" lt rgb "cyan", \
     'v1times_3' using 1:3 title "size3 Single RNG" lt rgb "green", \
     'v2times_3' using 1:3 title "size3 Distributed RNG" lt rgb "violet"


set term postscript eps enhanced color 
set title "Performance Comparison Sample Bucket Sort vs Naive Bucket Sort";
set output '|ps2pdf - Binningcomp.pdf';
set logscale y;
set xlabel "number of processes"; set ylabel "Total Completion Time (sec)";

plot 'v3times_1' using 1:($2+$3) title "size1 Sample Binning" lt rgb "red", \
     'v4times_1' using 1:($2+$3) title "size1 Naive Binning" lt rgb "black"


