/* 
 * Solves the Panfilov model using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory 
 * and reimplementation by Scott B. Baden, UCSD
 * 
 * Modified and  restructured by Didem Unat, Koc University
 * Modified for GPU by Thomas Chuna and Aaron Magilligan
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

// Timer: Make successive calls and take a difference to get the elapsed time.
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
void kernel(double* E, double* R, const int n, const double alpha, const double kk, 
		const double dt, const double a, const double epsilon, 
        	const double M1, const double M2, const double b)
{
	//kernel #1 to enforce boundary conditions

	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if(index <= n){
		E[ IMAT(index,0) ] = E[ IMAT(index,2) ];
		E[ IMAT(index,n+1) ] = E[ IMAT(index,n-1) ];
		E[ IMAT(0,index) ] = E[ IMAT(2,index) ];
		E[ IMAT(n+1,index) ] = E[ IMAT(n-1,index) ];
	}
	

	int i =  blockIdx.y*blockDim.y + threadIdx.y;
	int j =  blockIdx.x*blockDim.x + threadIdx.x;
	double tmpR =R[ IMAT(j,i) ]; double tmpR2 =0.0;
	double tmpE =0.0; double tmpE2 =0.0;

	__syncthreads();
	//kernel #2 to solve PDE
	if(i<=n && j <=n){
		tmpE = E[ IMAT(j,i) ]+alpha*(E[ IMAT(j,i+1) ]
			  + E[ IMAT(j,i-1) ]-4*E[ IMAT(j,i) ]
			  + E[ IMAT(j+1,i) ]+E[ IMAT(j-1,i) ]);
	}

	__syncthreads();
	//kernel #3 to solve E ODE
	if(i<=n && j <=n)
		tmpE2 = tmpE - dt*(kk*tmpE*( tmpE - a )*( tmpE - 1 ) + tmpE*tmpR );

	__syncthreads();
	//kernel #4 to solve R ODE
	if(i<=n && j <=n)
		tmpR2 = tmpR + dt*( epsilon + M1*tmpR / ( tmpE2 + M2) )*( -tmpR - kk*tmpE2*(tmpE2 - b - 1 ) );

	E[ IMAT(j,i) ] = tmpE2;
	R[ IMAT(j,i) ] = tmpR2;

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

	// Various constants - these definitions shouldn't change
	const double a = 0.1, b = 0.1, kk = 8.0, M1 = 0.07, M2 = 0.3, epsilon = 0.01, d = 5e-5;

	double T = 1000.0;
	int n = 200;
	int plot_freq = 0;
	int px = 1, py = 1;
	int no_comm = 0;
	int num_threads=1;

	cmdLine( argc, argv, T, n,px, py, plot_freq, no_comm, num_threads);
    int array_size = sizeof(double)*(n+2)*(n+2);

	// Allocate contiguous memory for solution arrays
	// The computational box is defined on [1:m+1,1:n+1]
	// We pad the arrays in order to facilitate differencing on the 
	// boundaries of the computation box

	double *E, *R;
	E = (double*)malloc(array_size);
	R = (double*)malloc(array_size);


	// Initialization
	for ( int j = 1; j <= n; j++ )
		for ( int i = 1; i <= n; i++ )
			E[ IMAT(j,i) ] = R[ IMAT(j,i) ] = 0;

	for ( int j = 1; j <= n; j++ )
	        for ( int i = n/2 + 1; i <= n; i++ )
			E[ IMAT(j,i) ] = 1.0;

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

	cout << "Grid Size       : " << n << endl; 
	cout << "Duration of Sim : " << T << endl; 
	cout << "Time step dt    : " << dt << endl; 
	cout << "Process geometry: " << px << " x " << py << endl;
	if (no_comm)
		cout << "Communication   : DISABLED" << endl;

	cout << endl;

	// Start the timer
	double t0 = getTime();


	// Simulated time is different from the integer timestep number

	// Simulated time
	double t = 0.0;

	// Integer timestep number
	int niter=0;
    
	//Initialize and Allocate GPU Memory
	double *d_E, *d_R;
	cudaMalloc((void**)&d_E, array_size);
	cudaMalloc((void**)&d_R, array_size);
	//Copy to CUDA kernels
	cudaMemcpy(d_E, E, array_size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_R, R, array_size, cudaMemcpyHostToDevice);

	//Define Grid and block sizes
	dim3 DimGrid(ceil((n+2)/16.0), ceil((n+2)/16.0), 1);
	dim3 DimBlock(16, 16, 1);

	while ( t < T ) {
	
		t += dt;
		niter++;
		kernel<<<DimGrid,DimBlock>>>(d_E, d_R,n,alpha,kk,dt,a,epsilon, M1, M2, b);

		if ( plot_freq ) {
			int k = ( (int) (t / plot_freq) );
			if ( (t - k*plot_freq ) < dt ) {
		                cudaMemcpy(E, d_E, array_size, cudaMemcpyDeviceToHost);
				splot(E,t,niter,n+2,n+2);
			}
		}
	} //end of while loop

	//Retrieve Result and Clean up GPU memory
	cudaMemcpy(E, d_E, array_size, cudaMemcpyDeviceToHost);
	cudaMemcpy(R, d_R, array_size, cudaMemcpyDeviceToHost);
	cudaFree(d_E);
	cudaFree(d_R);


	//Run Diagnostics
	double time_elapsed = getTime() - t0;
	double Gflops = (double)(niter * (1E-9 * n * n ) * 28.0) / time_elapsed;
	double BW = (double)(niter * 1E-9 * (n * n * sizeof(double) * 4.0  ))/time_elapsed;
	double mx;
	double l2norm = stats(E,n,&mx);

	//Output Diagnostics
	cout << "Number of Iterations        : " << niter << endl;
	cout << "Elapsed Time (sec)          : " << time_elapsed << endl;
	cout << "Sustained Gflops Rate       : " << Gflops << endl; 
	cout << "Sustained Bandwidth (GB/sec): " << BW << endl << endl; 
	cout << "Max: " << mx <<  " L2norm: "<< l2norm << endl;

	if (plot_freq)	{	
		cout << "\n\nEnter any input to close the program and the plot..." << endl;
		getchar();
	}

	//Clean Up CPU memory and close up shop
	free (E);
	free (R);
    
    
	return 0;
}
