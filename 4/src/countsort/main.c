#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>

#include "countsort_paralleli.h"
#include "countsort_parallelj.h"


int cmp_func(const void *a , const void *b)
{
    return ( *(int *) a - *(int *) b );
}


//
// Compares the arrays stored at a and b. Returns 1 if their elements are the same, 0 otherwise.
//
int difference( int * a, int * b, int n ) {

    int i;
    for(i = 0 ; i < n ; ++i){
        if(a[i] != b[i])
            return 0;
    }
    return 1;
}


//
// Main method. Runs tests and outputs results to stdout.
//
// Results are output on a single line with values separated by commas.
// Values are output in the following order:
//
//  1.  Number of threads.
//  2.  Runtime of sequential quicksort routine.
//  3.  Runtime of OpenMP i-parallel countsort routine.
//  4.  Value indicating correctness of OpenMP i-parallel countsort implementation.
//  5.  Runtime of OpenMP j-parallel countsort routine.
//  6.  Value indicating correctness of OpenMP j-parallel countsort implementation.
//
int main(int argc, char* argv[])
{

    int n, diff;
    int numberOfThreads;
    double time_i, time_j, time_quicksort;

    srand(time(NULL));

    // Get array size and number of threads from the command line
    if (argc > 1) {

        n = strtol(argv[1], NULL, 10);
        numberOfThreads = strtol(argv[2], NULL, 10);

    } else {

        fprintf( stderr, "\nWhen running this program, please include the array size and the number of threads on command line.\n");
        return 0;
    }

    int *a = malloc(n*sizeof(int));
    int *b = malloc(n*sizeof(int));
    int *c = malloc(n*sizeof(int));
    int *temp = malloc(n*sizeof(int));

    for ( int i = 0; i < n; ++i ) {

        c[i] = b[i] = a[i] = rand() % 1000;
        temp[i] = 0;
    }

    // Output number of threads.

    // Run quicksort library routine.
    {
        double start = omp_get_wtime();
        qsort(a, n, sizeof(int), cmp_func);
        double stop = omp_get_wtime();
        time_quicksort = stop - start;
    }

    // Run i-version of countsort.
    {
        double start = omp_get_wtime();
        countsort_paralleli( n, b, temp, numberOfThreads );
        double stop = omp_get_wtime();
        time_i = stop - start;
    }
    diff = difference(a,b,n);


    if ( diff == 0 )
        fprintf( stderr, "i-version of countsort produces incorrect result.\n" );

    // Run j-version of countsort.
    {
        double start = omp_get_wtime();
        countsort_parallelj( n, c, temp, numberOfThreads );
        double stop = omp_get_wtime();
        time_j = stop - start;
    }

    diff = difference(a,c,n);

    if ( diff == 0 )
        fprintf( stderr, "j-version of countsort produces incorrect result.\n" );
    double s_i = time_quicksort/time_i;
    double s_j = time_quicksort/time_j;
    printf("%d %d 1 %.4e %.4e %.4e %.4e \n", n, numberOfThreads,s_i, s_j, s_i/numberOfThreads, s_j/numberOfThreads );


    free(a);
    free(b);
    free(c);
    free(temp);
}
