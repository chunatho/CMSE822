#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>

#include "sequential.h"
#include "atomic.h"
#include "reduction.h"


// Error tolerance for checking correctness of approximation.
const double tol = 1e-2;

// Compute error in calculated value of pi.
double getDifference( double calculatedPi )
{
    return calculatedPi - 3.14159265358979323846;
}

//
// Main method. Runs tests and outputs results to stdout.
//
// Results are output on a single line with values separated by commas.
// Values are output in the following order:
//
//  1.  Number of Monte Carlo iterations.
//  2.  Number of OpenMP threads used for parallel calculations.
//  3.  Value of pi computed by sequential implementation.
//  4.  Error in value of pi computed by sequential implementation.
//  5.  Walltime in seconds required for sequential computation.
//  6.  Value of pi computed by OpenMP atomic implementation.
//  7.  Error in value of pi computed by OpenMP atomic implementation.
//  8.  Walltime in seconds required for OpenMP atomic computaton.
//  9.  Speedup of OpenMP atomic computation compared to sequential computation.
//  10. Value of pi computed by OpenMP reduction implementation.
//  11. Error in value of pi computed by OpenMP reduction implementation.
//  12. Walltime in seconds required for OpenMP reduction computaton.
//  13. Speedup of OpenMP reduction computation compared to sequential computation.
//
int main( int argc, char* argv[] )
{
    // Variable declarations
    double time_sequential, time_atomic, time_reduction;
    double pi_sequential = 0.0, pi_atomic = 0.0, pi_reduction = 0.0;

    int64_t iterations;
    int numberOfThreads;

    // Get number of iterations and number of threads from the command line
    if(argc > 1) {

        iterations = strtoimax(argv[1], NULL, 10);
        numberOfThreads = strtol(argv[2], NULL, 10);

    } else {

        fprintf( stderr, "\nWhen running this program, please include number of iterations and number of threads on command line.\n");
        return 0;
    }

    // Time sequential calculation
    {
        double start = omp_get_wtime();
        pi_sequential = sequentialCompute(iterations);
        double stop = omp_get_wtime();
        time_sequential = stop - start;

    }

    // Time parallel calculation with atomics
    {
        double start = omp_get_wtime();
        pi_atomic = parallelCompute_atomic(iterations, numberOfThreads);
        double stop = omp_get_wtime();
        time_atomic = stop - start;
    }

    // Time parallel calculation with reduction
    {
        double start = omp_get_wtime();
        pi_reduction = parallelCompute_reduction(iterations, numberOfThreads);
        double stop = omp_get_wtime();
        time_reduction = stop - start;
    }

    // Compute errors in approximations of pi.
    double error_sequential = getDifference(pi_sequential);
    double error_atomic = getDifference(pi_atomic);
    double error_reduction = getDifference(pi_reduction);

//    fprintf(stderr, "sequential_pi is %E error  is: %E \n", pi_sequential, error_sequential );

//    fprintf(stderr, "atomic_pi is %E  error is : %E \n", pi_atomic,error_atomic );

//    fprintf(stderr, "reduction_pi is %E error  is : %E \n", pi_reduction,error_reduction );

    // Print iteration and thread counts.
   double s_r = time_sequential/time_reduction;
   double s_a = time_sequential/time_atomic;


    printf("%" PRId64 " %d 1 %.4e %.4e %.4e %.4e \n", iterations, numberOfThreads, s_r,s_a, s_r/numberOfThreads, s_a/numberOfThreads);
    
    if ( fabs( error_sequential ) > tol )
        fprintf( stderr, "%" PRId64 " Sequential calculation is INVALID -  error is = %E \n", iterations, error_sequential);

    if ( fabs( error_atomic ) > tol )
        fprintf( stderr, "%" PRId64 " Parallel atomic calculation is INVALID - error is = %E \n", iterations, error_atomic );

    if ( fabs( error_reduction ) > tol )
        fprintf( stderr, "%" PRId64 " Parallel atomic calculation is INVALID - error is = %E \n", iterations, error_reduction);

    return 0;
}

