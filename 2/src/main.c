//
//  Created by H. Metin Aktulga on 2/2/17.
//  Modified by Michael M. Crockatt on 9/8/17.
//  Copyright © 2017 Michigan State University. All rights reserved.
//

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <sys/time.h> // for clock_gettime()

#include "naiveTranspose.h"
#include "optTranspose.h"

// Method declarations
double getElapsed(struct timeval*, struct timeval*);

void initMatrix(double*, int, int);
int checkDifference(double*, double*, int, int);

// Main method
int main(int argc, char* argv[])
{
    // Variable declarations
    double naiveTime, optTime;
    double* naiveinput, *naiveoutput, *optinput, *optoutput;
    int  N, M, B;
    int  testing = 0, ninps = 0, curtest = 0;
    FILE *fin;
    struct timeval start, end;

    // Get the problem size and blocking factor from the command line
    if (argc == 3 || argc == 5) {

        if (!strcmp(argv[1],"test")) {
            testing = 1;
            fin = fopen(argv[2], "r");
            fscanf(fin, " %d", &ninps);
        } else if (!strcmp(argv[1],"perf")) {
            testing = 0;
        } else {
            fprintf(stderr, "\nERROR: Unknown option %s\n", argv[1]);
        }

    } else {
        fprintf(stderr, "\nERROR: This program can be executed in the \
\"test\" or \"perf\" modes:\n\n");
        fprintf(stderr,"           transpose test inputfile\n");
        fprintf(stderr,"           transpose perf N M B\n\n");
        fprintf(stderr,"           inputfile format: starts with a line indicating \
the number of tests,\n");
        fprintf(stderr,"                             followed the tests specified with \
N M B values.\n\n");
        fprintf(stderr,"           N: number of matrix rows\n");
        fprintf(stderr,"           M: number of matrix columns\n");
        fprintf(stderr,"           B: blocking factor.\n");
        return 0;
    }

    // If doing testing run, print header for output to CSV file.
    if (testing) printf( "N,M,B,naive_time,opt_time,opt_correct\n" );

    do {
        // main loop: executed only once for performance testing to ensure cold start.
        //            However, it is okay to do several accuracy checks (test runs)
        //            at once in a single execution.

        if (testing) {
            curtest++;
            // read from the input file
            fscanf(fin, " %d %d %d", &N, &M, &B);
        } else {
            // otherwise get the command line arguments
            N = strtol(argv[2], NULL, 10);
            M = strtol(argv[3], NULL, 10);
            B = strtol(argv[4], NULL, 10);
        }

        if (testing) {
            printf( "%d,%d,%d,", N, M, B );
        } else {
            printf("ninps = %d, curtest = %d\n", ninps, curtest);
            printf("N = %d, M = %d, B = %d\n", N, M, B);
        }

        // Allocations
        if (testing) {
            naiveinput = (double *) malloc(N*M*sizeof(double));
            initMatrix(naiveinput, N, M);
            optinput = naiveinput;
        } else {
            naiveinput = (double *) malloc(N*M*sizeof(double));
            optinput = (double *) malloc(N*M*sizeof(double));
        }

        naiveoutput = (double *) malloc(M*N*sizeof(double));
        optoutput = (double *) malloc(M*N*sizeof(double));

        // Time naive calculation
        gettimeofday(&start, NULL);
        naiveTranspose(naiveoutput, naiveinput, N, M);
        gettimeofday(&end, NULL);
        naiveTime = getElapsed(&start, &end);

        // Time opt calculation with atomics
        gettimeofday(&start, NULL);
        optTranspose(optoutput, optinput, N, M, B);
        gettimeofday(&end, NULL);
        optTime = getElapsed(&start, &end);

        // How do results compare?
        if (testing) {
            printf( "%.4e,%.4e,%d\n", naiveTime, optTime,
                    checkDifference(naiveoutput, optoutput, M, N) );
        } else {
            printf("Naive time: %.3f sec\n", naiveTime);
            printf("Opt. time:  %.3f sec\n", optTime);
            printf("Achieved speedup= %.2f\n\n", naiveTime/optTime);
        }

        if (naiveinput) free(naiveinput);
        if (naiveoutput) free(naiveoutput);
        if (optinput && !testing) free(optinput);
        if (optoutput) free(optoutput);

    } while(curtest < ninps);

    if (testing) fclose(fin);

    return 0;
}


double getElapsed(struct timeval *start, struct timeval *end)
{
    double secs_used=(end->tv_sec - start->tv_sec); //avoid overflow by subtracting first
    double micros_used= ((secs_used*1000000) + end->tv_usec) - (start->tv_usec);

    return micros_used/1000000.0;
}


void initMatrix(double *matrix, int N, int M)
{
    srand(time(NULL));

    for ( int i = 0; i < N; ++i ) {
        for ( int j = 0; j < M; ++j ) {
            matrix[i*M+j] = rand()*1.0/INT_MAX;
        }
    }
}


int checkDifference(double *matrix1, double *matrix2, int N, int M)
{
    for ( int i = 0; i < N; ++i ) {
        for ( int j = 0; j < M; ++j ) {
            if (matrix1[i*M+j] != matrix2[i*M+j]) {
                fprintf(stderr, "matrices differ at position [%d,%d]: ", i, j);
                fprintf(stderr, "matrix1[%d,%d] = %.5f, ", i, j, matrix1[i*M+j]);
                fprintf(stderr, "matrix2[%d,%d] = %.5f\n", i, j, matrix2[i*M+j]);
                return 0;
            }
        }
    }

    return 1;
}
