
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <stdio.h>

//
// TODO: - Parallelize the *** i-loop *** in this source file.
//
//everybody gets a cop of a then you break a into pieces and have each thread run its counting
//schemes for the elementes its been given. you concatinate the increment lists
int64_t count;
void countsort_paralleli( const int n, int * const a, int * const temp, const int numberOfThreads ) {
   int thread_count = numberOfThreads;
#pragma omp parallel for private(count) num_threads(thread_count)
      for ( int i = 0 ; i < n ; i++ ) {
         count = 0;
         for ( int j = 0 ; j < n ; j++ ) {
            if(a[j] < a[i])
               count++;
            else if((a[j] == a[i]) && (j < i))
               count++;
      }
//#pragma omp critical
      temp[count] = a[i];

//      memcpy(a[start], temp[start], length*sizeof(int) );
   }
#pragma omp barrier
//   printf("starting paralleli memcpy \n");

   memcpy( a, temp, n * sizeof(int) );
}
