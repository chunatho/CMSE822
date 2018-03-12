#!/usr/bin/gnuplot

set term postscript eps enhanced color 
set output '|ps2pdf - cs100comparison.pdf';
set title "Comparison of 10^2 countsort loop parallelization";
set xlabel "number of threads"; set ylabel "absolute speedup";
set logscale y;
plot 'cs1' using 2:3 title "Sequential" lt rgb "red", \
     'cs1' using 2:4 title "paralleli" lt rgb "black", \
     'cs1' using 2:5 title "parallelj" lt rgb "blue"

set output '|ps2pdf - cs1000comparison.pdf';
set title "Comparison of 10^3 countsort loop parallelization";
set xlabel "number of threads"; set ylabel "absolute speedup";
set logscale y;
plot 'cs2' using 2:3 title "Sequential" lt rgb "red", \
     'cs2' using 2:4 title "paralleli" lt rgb "black", \
     'cs2' using 2:5 title "parallelj" lt rgb "blue"

set output '|ps2pdf - cs10000comparison.pdf';
set title "Comparison of 10^4 countsort loop parallelization";
set xlabel "number of threads"; set ylabel "absolute speedup";
set logscale y;
plot 'cs3' using 2:3 title "Sequential" lt rgb "red", \
     'cs3' using 2:4 title "paralleli" lt rgb "black", \
     'cs3' using 2:5 title "parallelj" lt rgb "blue"




set output '|ps2pdf - pi10-5comparison.pdf';
set title "Comparison of 10^5 speed up across parallelization";
set xlabel "thread number"; set ylabel "relative speed up";

plot 'pi1' using 2:3 title "sequential", \
     'pi1' using 2:4 title "reduction", \
     'pi1' using 2:5 title "atomic"
set output '|ps2pdf - pi10-7comparison.pdf';
set title "Comparison of 10^7 speed up across parallelization";
set xlabel "thread number"; set ylabel "relative speed up";

plot 'pi2' using 2:3 title "sequential", \
     'pi2' using 2:4 title "reduction", \
     'pi2' using 2:5 title "atomic"
set output '|ps2pdf - pi10-9comparison.pdf';
set title "Comparison of 10^9 speed up across parallelization";
set xlabel "thread number"; set ylabel "relative speedup";

plot 'pi3' using 2:3 title "sequential", \
     'pi3' using 2:4 title "reduction", \
     'pi3' using 2:5 title "atomic"


set output '|ps2pdf - pi10-5effcomp.pdf';
set title "Comparison of 10^5 parallelization efficency";
set xlabel "number of threads"; set ylabel "Efficiency";
plot 'pi1' using 2:6 title "reduction" lt rgb "black", \
     'pi1' using 2:7 title "atomic" lt rgb "blue"

set output '|ps2pdf - pi10-7effcomp.pdf';
set title "Comparison of 10^7 parallelization efficency";
set xlabel "number of threads"; set ylabel "Efficency";
plot 'pi2' using 2:6 title "reduction" lt rgb "black", \
     'pi2' using 2:7 title "atomic" lt rgb "blue"

set output '|ps2pdf - pi10-9effcomp.pdf';
set title "Comparison of 10^9 parallelization efficency";
set xlabel "number of threads"; set ylabel "Efficiency";
plot 'pi3' using 2:6 title "reduction" lt rgb "black", \
     'pi3' using 2:7 title "atomic" lt rgb "blue"

