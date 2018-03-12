#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
//#include <random>

//
// TODO: You need to implement a sequential estimation for PI using the Monte Carlo method.
//       Use rand_r() here as well for fair comparison to OpenMP parallel versions.
//
double sequentialCompute( int64_t iterations )
{

int64_t Ntoss = iterations;
int64_t Nhit = 0;
double x, y, r2, pi;

int seed = 17;
srand(seed);

for (int toss = 0; toss < Ntoss; toss++) {
//   x = Random::get(0.f,1.f)
   x = drand48();
   y = drand48();
   r2 = x * x + y * y;
   if (r2 <= 1) {Nhit++;}
}
   pi = 4 * Nhit/((double) Ntoss);
   return pi;
}
