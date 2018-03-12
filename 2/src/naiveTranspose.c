//
//  Created by H. Metin Aktulga on 2/2/17.
//  Modified by Michael M. Crockatt on 9/8/17.
//  Copyright © 2017 Michigan State University. All rights reserved.
//

#include "naiveTranspose.h"

// Naive transpose operation
void naiveTranspose(double *output, double *input, int N, int M)
{
    for ( int i = 0; i < N; ++i ) {
        for ( int j = 0; j < M; ++j ) {
            output[j*N+i] = input[i*M+j];
        }
    }
}
