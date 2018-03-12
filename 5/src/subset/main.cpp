
//============================================================================//
//============================================================================//
//============================================================================//
//                                                                            //
//                DO NOT MODIFY THE CONTENTS OF THIS FILE.                    //
//                                                                            //
//============================================================================//
//============================================================================//
//============================================================================//


#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <omp.h>
#include <time.h>

#include "generate_subset.hpp"


// Holds the current minimum value found so far.
int global_min = std::numeric_limits<int>::max();

// Sum of all of the weights.
int total_sum;


//------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------
int main( int argc, char* argv[] ) {

    double start_time, end_time;
    int nweights, nthreads;
    std::vector<int> ans, S;

    // Get thread and weight counts from command line.
    if ( argc == 3 ) {

        nweights = strtol(argv[1], NULL, 10);
        nthreads = strtol(argv[2], NULL, 10);

    } else {

        printf( "\nThis program requires two input arguments: "
                "the number of weights and the number of threads.\n" );
        return 0;
    }

    srand(12);

    for ( int i = 0; i < nweights; i++ ) {

        S.push_back( rand() % 100 );
        total_sum += S[i];
    }

    printf("Weight set: ");

    for ( int i = 0; i < S.size(); i++)
        printf("%d ", S[i]);

    printf("Total sum of weights: %d \n", total_sum);

    start_time = omp_get_wtime();

    #pragma omp parallel num_threads(nthreads) shared(S)
    {
        #pragma omp single nowait
        GenerateSubset(0, ans, S);
    }

    end_time = omp_get_wtime();

//    printf("\nMinimum difference found: %d\n", global_min);
    printf("Time %d %.4e\n",nthreads,  end_time - start_time);
}
