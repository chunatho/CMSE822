#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

//
// TODO: You need to implement an OpenMP parallel version using the OpenMP
//       reduction clause.
//       Use rand_r() for thread safe random number generation.
//       More details about rand_r() is here: http://linux.die.net/man/3/rand_r
//       You must also make sure that each thread start with a DIFFERENT SEED!
//


double parallelCompute_reduction( int64_t iterations, int numberOfThreads )
{
int thread_count = numberOfThreads;
int64_t global_Nhit = 0;
double pi;
unsigned int seed[thread_count];

for(int i =0; i < thread_count; i++) {
   seed[i]=i;
}
#pragma omp parallel reduction(+: global_Nhit) num_threads(thread_count)
{
int thread_num = omp_get_thread_num();


#pragma omp for

   for(int64_t i = 0; i < iterations;i++){
      double x = (double)rand_r(&seed[thread_num])/(RAND_MAX);
      double y = (double)rand_r(&seed[thread_num])/(RAND_MAX);
      double r2 = x * x + y * y;
      if (r2 <= 1)
         {global_Nhit++;}
   }

//#pragma omp barrier
}
pi = 4* global_Nhit/((double) iterations);
return pi;

}
