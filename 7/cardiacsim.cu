/* 
 * Solves the Panfilov model using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory 
 * and reimplementation by Scott B. Baden, UCSD
 * 
 * Modified and  restructured by Didem Unat, Koc University
 *
 */
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <math.h>
#include <sys/time.h>
using namespace std;

# define IMAT(j,i) ( (i) + (n+2)*(j) ) 

// Utilities
// 

// Timer
// Make successive calls and take a difference to get the elapsed time.
static const double kMicro = 1.0e-6;
double getTime() {
	struct timeval TV;
    struct timezone TZ;

    const int RC = gettimeofday(&TV, &TZ);
    if (RC == -1) {
            cerr << "ERROR: Bad call to gettimeofday" << endl;
            return(-1);
    }

    return ( ( (double) TV.tv_sec ) + kMicro*( (double) TV.tv_usec) );

}  // end getTime()

// Allocate a 2D array
double ** alloc2D(int n) {

	// Comments here would be nice
	double **E;

	E = (double**)malloc(sizeof(double*)*n + sizeof(double)*n*n);
	assert(E);

	for( int j = 0; j < n; j++ ) 
		E[j] = (double*)(E+n) + j*n;
	return(E);
}
    
// Reports statistics about the computation
// These values should not vary (except to within roundoff)
// when we use different numbers of  processes to solve the problem
double stats(double *E, int n, double *_mx) {
	double mx = -1;
	double l2norm = 0;

	for ( int j = 1; j <= n; j++ ) 
		for ( int i = 1; i <= n; i++ ) {

			l2norm += E[ IMAT(j,i) ]*E[ IMAT(j,i) ];
			if ( E[ IMAT(j,i) ] > mx )
				mx = E[ IMAT(j,i) ];

		}	

	*_mx = mx;
	l2norm /= ( (double) (n*n) );
	l2norm = sqrt(l2norm);

	return l2norm;
 }

// External functions
extern "C" {
    void splot(double *E, double T, int niter, int m, int n);
}
void cmdLine(int argc, char *argv[], double& T, int& n, int& px, int& py, int& plot_freq, int& no_comm, int&num_threads);


__global__
void first_kernel(double* E_prev, const int n)
{
    int j = blockDim.x * blockIdx.x + threadIdx.x;
    if(j<=n) {
        E_prev[ IMAT(j,0)  ] = E_prev[ IMAT (j,2) ];

		E_prev[ IMAT(j,n+1) ] = E_prev[ IMAT(j,n-1) ];

		E_prev[ IMAT(0,j) ] = E_prev[ IMAT(2,j) ];

		E_prev[ IMAT(n+1,j) ] = E_prev[ IMAT(n-1,j) ];
    }
}
__global__
void second_kernel(double* E, double* E_prev, const int n, const double alpha)
{
    int i = threadIdx.x;
    int j = blockIdx.x;
    int k = blockDim.x * blockIdx.x + threadIdx.x;
    if(k<=(n+2)*(n+2)) {
			E[ IMAT(j,i) ] = E_prev[ IMAT(j,i) ]+alpha*(E_prev[ IMAT(j,i+1) ]
					  + E_prev[ IMAT(j,i-1) ]-4*E_prev[ IMAT(j,i) ]
					  + E_prev[ IMAT(j+1,i) ]+E_prev[ IMAT(j-1,i) ]);
    }
}

__global__
void third_kernel(double* E, double* R, const int n, const double kk, const double dt, const double a)
{
    int i = threadIdx.x;
    int j = blockIdx.x;
    int k = blockDim.x * blockIdx.x + threadIdx.x;
    if(k<=(n+2)*(n+2)) {
			E[ IMAT(j,i) ] = E[ IMAT(j,i) ] - dt*(kk*E[ IMAT(j,i) ]*( E[ IMAT(j,i) ] - a )*( E[ IMAT(j,i) ] - 1 ) 
                + E[ IMAT(j,i) ]*R[ IMAT(j,i) ] );
    }
}
__global__
void fourth_kernel(double* E, double* R, const int n, const double kk, const double dt, const double epsilon, 
        const double M1, const double M2, const double b)
{
    int i = threadIdx.x;
    int j = blockIdx.x;
    int k = blockDim.x * blockIdx.x + threadIdx.x;
    if(k<=(n+2)*(n+2)) {
			R[ IMAT(j,i) ] = R[ IMAT(j,i) ] + dt*( epsilon + M1*R[ IMAT(j,i) ] / ( E[ IMAT(j,i) ] + M2) )*( -R[ IMAT(j,i) ]
                 - kk*E[ IMAT(j,i) ]*(E[ IMAT(j,i) ] - b - 1 ) );
    }
}
void simulate (double* d_E,  double* d_E_prev,double* d_R,
   			   const double alpha, const int n, const double kk,
   			   const double dt, const double a, const double epsilon,
               const double M1,const double  M2, const double b) {
/*
	double * E = (double*)(E_old + n + 2);
	double * E_prev = (double*)(E_prev_old + n + 2);
	double * R = (double*)(R_old + n + 2);
*/
	/* 
	* Copy data from boundary of the computational box 
	* to the padding region, set up for differencing
	* on the boundary of the computational box
	* Using mirror boundaries
	*/
//    printf("ceil((n+2)/256.0): %f \t n: %d \n", ceil((n+2)/256.0), n);
    
    first_kernel<<<ceil((n+2)),n+2>>>(d_E_prev,n);
    second_kernel<<<ceil(n+2),n+2>>>(d_E,d_E_prev,n,alpha);
    third_kernel<<<ceil(n+2),n+2>>>(d_E,d_R,n,kk,dt,a);
    fourth_kernel<<<ceil(n+2),n+2>>>(d_E,d_R,n,kk,dt,epsilon,M1,M2,b);
    

/*
	for ( int j = 1; j <= n; j++ ){ 
		E_prev[ IMAT(j,0) ] = E_prev[ IMAT(j,2) ];
        if(j==27) printf("j: %d \t h_E_prev: %f \t E_prev: %f \n", j, h_E_prev[j], E_prev[j]);
    }
	for ( int j = 1; j <= n; j++ ) 
		E_prev[ IMAT(j,n+1) ] = E_prev[ IMAT(j,n-1) ];

	for ( int i = 1; i <= n; i++ ) 
		E_prev[ IMAT(0,i) ] = E_prev[ IMAT(2,i) ];

	for ( int i = 1; i <=n ; i++ ) 
		E_prev[ IMAT(n+1,i) ] = E_prev[ IMAT(n-1,i) ];
*/
/*
	// Solve for the excitation, the PDE
	for ( int j = 1; j <= n; j++ ){
		for ( int i = 1; i <= n; i++ ) {
			E[ IMAT(j,i) ] = E_prev[ IMAT(j,i) ]+alpha*(E_prev[ IMAT(j,i+1) ]
					  + E_prev[ IMAT(j,i-1) ]-4*E_prev[ IMAT(j,i) ]
					  + E_prev[ IMAT(j+1,i) ]+E_prev[ IMAT(j-1,i) ]);
		}	
}
*/
	/* 
	* Solve the ODE, advancing excitation and recovery to the
	* next timtestep
	*/
/*
	for ( int j = 1; j <= n; j++ ) {
		for ( int i=1; i<=n; i++)
			E[ IMAT(j,i) ] = E[ IMAT(j,i) ] - dt*(kk*E[ IMAT(j,i) ]*( E[ IMAT(j,i) ] - a )*( E[ IMAT(j,i) ] - 1 ) + E[ IMAT(j,i) ]*R[ IMAT(j,i) ] );
	}
*/
/*
	for ( int j = 1; j <= n; j++ ) {
		for ( int i = 1; i <= n; i++ )
			R[ IMAT(j,i) ] = R[ IMAT(j,i) ] + dt*( epsilon + M1*R[ IMAT(j,i) ] / ( E[ IMAT(j,i) ] + M2) )*( -R[ IMAT(j,i) ] - kk*E[ IMAT(j,i) ]*(E[ IMAT(j,i) ] - b - 1 ) );
	}	
*/
}

// Main program
int main (int argc, char** argv) {

	/*
	*  Solution arrays
	*   E is the "Excitation" variable, a voltage
	*   R is the "Recovery" variable
	*   E_prev_old is the Excitation variable for the previous timestep,
	*      and is used in time integration
	*/
//	double **E_old, **R_old, **E_prev_old;

	// Various constants - these definitions shouldn't change
	const double a = 0.1, b = 0.1, kk = 8.0, M1 = 0.07, M2 = 0.3, epsilon = 0.01, d = 5e-5;

	double T = 1000.0;
	int n = 200;
	int plot_freq = 0;
	int px = 1, py = 1;
	int no_comm = 0;
	int num_threads=1; 

	cmdLine( argc, argv, T, n,px, py, plot_freq, no_comm, num_threads);

	// Allocate contiguous memory for solution arrays
	// The computational box is defined on [1:m+1,1:n+1]
	// We pad the arrays in order to facilitate differencing on the 
	// boundaries of the computation box
//	E_old = alloc2D( n + 2 );
//	E_prev_old = alloc2D( n + 2 );
//	R_old = alloc2D( n + 2 );

    double *E; //= (double*)(E_old + n + 2);
    double *E_prev; // = (double*)(E_prev_old + n + 2);
    double *R; // = (double*)(R_old + n + 2);
	
    E = (double *)malloc((n+2)*(n+2)*sizeof(double));
    E_prev = (double *)malloc((n+2)*(n+2)*sizeof(double));
    R = (double *)malloc((n+2)*(n+2)*sizeof(double));

    // Initialization
	for ( int j = 1; j <= n; j++ )
		for ( int i = 1; i <= n; i++ )
			E_prev[ IMAT(j,i) ] = R[ IMAT(j,i) ] = 0;

	for ( int j = 1; j <= n; j++ )
		for ( int i = n/2 + 1; i <= n; i++ )
			E_prev[ IMAT(j,i) ] = 1.0;

	for ( int j = n/2 + 1; j <= n; j++ )
		for ( int i = 1; i <= n; i++ )
			R[ IMAT(j,i) ] = 1.0;

	double dx = 1.0/n;

	// For time integration, these values shouldn't change 
	double rp= kk*( b + 1 )*( b + 1 ) / 4;
	double dte = ( dx*dx ) / ( d*4 +  ( dx*dx ) * ( rp + kk ) );
	double dtr = 1 / ( epsilon + (M1/M2)*rp );
	double dt = ( dte < dtr ) ? 0.95*dte : 0.95* dtr;
	double alpha = d*dt / ( dx*dx );
/*
	cout << "Grid Size       : " << n << endl; 
	cout << "Duration of Sim : " << T << endl; 
	cout << "Time step dt    : " << dt << endl; 
	cout << "Process geometry: " << px << " x " << py << endl;
	if (no_comm)
		cout << "Communication   : DISABLED" << endl;

	cout << endl;
*/
	// Start the timer
	double t0 = getTime();


	// Simulated time is different from the integer timestep number

	// Simulated time
	double t = 0.0;

	// Integer timestep number
	int niter=0;
    double *d_E_prev, *d_E, *d_R;
    cudaMalloc(&d_E,(n+2)*(n+2)*sizeof(double));
    cudaMalloc(&d_E_prev,(n+2)*(n+2)*sizeof(double));
    cudaMalloc(&d_R,(n+2)*(n+2)*sizeof(double));
    cudaMemcpy(d_E, E, (n+2)*(n+2)*sizeof(double), cudaMemcpyHostToDevice); 
    cudaMemcpy(d_E_prev, E_prev, (n+2)*(n+2)*sizeof(double), cudaMemcpyHostToDevice); 
    cudaMemcpy(d_R, R, (n+2)*(n+2)*sizeof(double), cudaMemcpyHostToDevice); 

	while ( t < T ) {
	
		// Don't know what the fuck this is supposed to do since there's no comment
		t += dt;
		niter++;


		simulate(d_E, d_E_prev, d_R, alpha, n, kk, dt, a, epsilon, M1, M2, b); 

		//swap current E with previous E
		double *tmp = d_E; d_E = d_E_prev; d_E_prev = tmp;

		if ( plot_freq ) {
			int k = ( (int) (t / plot_freq) );
			if ( (t - k*plot_freq ) < dt ) {
				splot(d_E,t,niter,n+2,n+2);
			}
		}
	} //end of while loop
    cudaMemcpy(E, d_E, (n+2)*(n+2)*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(E_prev, d_E_prev, (n+2)*(n+2)*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(R, d_R, (n+2)*(n+2)*sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_E);
    cudaFree(d_E_prev);
    cudaFree(d_R);

	double time_elapsed = getTime() - t0;

	double Gflops = (double)(niter * (1E-9 * n * n ) * 28.0) / time_elapsed;
	double BW = (double)(niter * 1E-9 * (n * n * sizeof(double) * 4.0  ))/time_elapsed;
/*
	cout << "Number of Iterations        : " << niter << endl;
	cout << "Elapsed Time (sec)          : " << time_elapsed << endl;
	cout << "Sustained Gflops Rate       : " << Gflops << endl; 
	cout << "Sustained Bandwidth (GB/sec): " << BW << endl << endl; 
*/
	double mx;
	double l2norm = stats(E_prev,n,&mx);
//	cout << "Max: " << mx <<  " L2norm: "<< l2norm << endl;
    printf("v1 \t %d \t %.10e \t %.10e \t %.6e \n",n,mx,l2norm,time_elapsed);
    if (plot_freq)	{	
		cout << "\n\nEnter any input to close the program and the plot..." << endl;
		getchar();
	}

/*
	free (E_old);
	free (E_prev_old);
	free (R_old);
*/ 
       
	free (E);
	free (E_prev);
	free (R);
	
    return 0;
}
