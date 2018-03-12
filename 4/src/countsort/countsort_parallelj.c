#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <stdio.h>

//
// TODO: - Parallelize the *** j-loop *** in this source file.
//
//for j parallel you parallelize the comparisons so that given an element each thread
//runs comparisons on subset of the list, you agregae incrementation at the end

void counter2(int start, int stop,int i,int* global_count, int * const a);

void countsort_parallelj( const int n, int * const a, int * const temp, const int numberOfThreads ) {
   int thread_count = numberOfThreads;
   int start, stop;
   int length = (n-n%thread_count+thread_count)/thread_count;
   int global_count=0;

//   printf("starting parallelj i for loop \n");
   for ( int i = 0 ; i < n ; i++ ) {
      global_count = 0;

#pragma omp parallel for reduction(+:global_count) num_threads(thread_count)
      for ( int k = 0; k < thread_count; k++) {
         start = length*k;
         if (length*(k+1) < n){stop = length*(k+1);}
         else {stop = n;}
//         printf(" start and stop %d %d \n", start, stop);
         for ( int j = start ; j < stop ; j++ ) {
         if(a[j] < a[i])
            global_count++;
         else if((a[j] == a[i]) && (j < i))
            global_count++;
         }
      }
#pragma omp barrier
//      printf("global_count %d \n",global_count);
      temp[global_count] = a[i];
   }
   // TODO: Modify the code below so that the copy can be made OpenMP parallel
//#pragma omp parallel num_threads(thread_count)
//   for(int k = 0; k < thread_count; k++){
//      start = length*k;
//      if (length*(k+1) < n){stop = length*(k+1);}
//      else {stop = n;} 
//      memcpy(a[start], temp[start], length * sizeof(int) );
//   }

//   printf("starting parallelj memcpy \n");
   memcpy( a, temp, n * sizeof(int) );
}
