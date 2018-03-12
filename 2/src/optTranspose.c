//
//  Created by H. Metin Aktulga on 2/2/17.
//  Modified by Michael M. Crockatt on 9/8/17.
//  Copyright © 2017 Michigan State University. All rights reserved.
//

# include <math.h>

#include "optTranspose.h"

//
// TODO: You need to implement this function.
//
// Cache blocked transpose operation
//
// If something doesn't work switch the N and the M in the outpt line
// and double check that *output / *input is the right way 
// to reference those arguements

void optTranspose(double *output, double *input, int N, int M, int B)
{
 int kMax, lMax;
 for (int i = 0; i < N; i += B) {
    for (int j = 0; j < M; j += B) {
       // transpose the block beginning at [i,j]
       kMax =fmin(N,i+B);
       lMax =fmin(M,j+B);
       for (int k = i; k < kMax; ++k) {
          for (int l = j; l < lMax; ++l) {
               output[k + l*N] = input[l + k*M];
          }
       }
    }
 } 
}
