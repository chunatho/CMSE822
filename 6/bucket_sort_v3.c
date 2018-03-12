#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

// Comparison function used by qsort
int compare_dbls(const void* arg1, const void* arg2)
{
    double a1 = *(double *) arg1;
    double a2 = *(double *) arg2;
    if (a1 < a2) return -1;
    else if (a1 == a2) return 0;
    else return 1;
}
// Sort the array in place
void qsort_dbls(double *array, int array_len)
{
    qsort(array, (size_t)array_len, sizeof(double), compare_dbls);
}

void main(int argc, char* argv[]){
//Launch Message Passing Interface
    MPI_Init(NULL,NULL);
    int p, my_rank,S;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    double *input_array;   // the input array
    double *bucketlist;    //this array will contain input elements in order of the processors
                            //e.g elements of process 0 will be stored first, then elements of process 1, and so on
    double *local_array;    //This array will contain the elements in each process
    double *orderedlist;
    double *S_local;
    double *S_sorted;
    double tstart1, tend1, tstart2, tend2;
    int *local_size_array;

    int n,i;

    int *scounts;           //This array will contain the counts of elements each processor will receive
    int *Sscounts;
    int *dspls;             //The relative offsets in bucketlist array where the elements of different processes
                            //will be stored
    int *local_dspls;
    int *S_dspls;
    int *bin_elements;      //it will keep track of how many elements have been included in the pth bi 
    n = atoi(argv[1]);
    S = (int) 12*p*log(n);
    //printf("S: %d \n",S);

    int stop = n/p;
    int Sstop = S/p;
    //Redefine S here to avoid the modulus problem, it's arbitrary anyway
    S = Sstop*p;

//Allocate memory
    tstart1= MPI_Wtime();
    if( my_rank == 0){
        stop += n%p;
        }
    input_array = malloc(stop*sizeof(double));
    orderedlist = malloc(n*sizeof(double));
    bucketlist = malloc(stop*sizeof(double));
    S_local = malloc(Sstop*sizeof(double));
    S_sorted = malloc(S*sizeof(double));
    scounts = malloc(p*sizeof(int));
    Sscounts = malloc(p*sizeof(int));
    local_size_array = malloc(p*sizeof(int));
    dspls = malloc(p*sizeof(int));
    local_dspls = malloc(p*sizeof(int));
    S_dspls = malloc(p*sizeof(int));
    bin_elements = malloc(p*sizeof(int));

// Every process generates random data
    for(i = 0 ; i < stop ; i++){
        double tmp = (double) rand() / RAND_MAX;
        input_array[i] = tmp * tmp;
        //printf("element %d :   %f     %f\n",i,input_array[i],tmp);
    }
//select S/p elements from input_array at random
    int randspot;

    for(i=0 ; i < Sstop ; i++){
        randspot = rand()%Sstop;
        S_local[i] = input_array[i];
    }
    S_dspls[0] = 0;
    for(i = 0 ; i< p-1 ;i++){
        S_dspls[i+1] = S_dspls[i] + Sstop;
    }

    for(i = 0 ; i < p ; i++){
        Sscounts[i] = Sstop ;
    }


    MPI_Gatherv(S_local, Sstop, MPI_DOUBLE, S_sorted, Sscounts, S_dspls, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double *pivots;
    pivots = malloc((p)*sizeof(double));
    pivots[p-1] = 1.0;
    if(my_rank ==0){
      for( i = 0 ; i < S ; i++){
      //printf("element %d: %f \n",i,S_sorted[i]);
      }
      qsort_dbls(S_sorted, S);
      for( i = 0 ; i < S ; i++){
      //printf("element %d: %f \n",i,S_sorted[i]);
      }
      for( i = 1 ; i < p ; i++){
      int tmp = i*S/p;
      printf("break %d at element %d: %f \n",i,tmp,S_sorted[tmp]);
      pivots[i-1] = S_sorted[tmp];
      }
    }
    
    MPI_Bcast(pivots,p, MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);

/*
    for( i = 0 ; i < p ; i++){
      printf("rank %d: pivot %d: %f \n",my_rank,i,pivots[i]);
      }
*/

    for(i = 0 ; i < p ; i++){
        scounts[i] = 0 ;
    }

//counting the elements in each processor
    double condition;
    for(i = 0 ; i < stop ; i++){
        int j = 0;
        condition = input_array[i];
        while(condition > 0){
          condition = input_array[i] -  pivots[j];
          j++;
        }
        scounts[j-1]++;
    }
    //Make sure everyone has the right counts using the above sorting
    /*
    for(i = 0 ; i < p ; i++){
        printf("scounts[%d] = %d\n",i,scounts[i]) ;
    }
    */


    for(i = 0 ; i<p ; i++){
        bin_elements[i] = scounts[i];
    }


    dspls[0] = 0;
    for(i = 0 ; i< p-1 ;i++){
        dspls[i+1] = dspls[i] + scounts[i];
    }

    int bin;
    int pos;
    for(i = 0 ; i < stop ; i++){
        int j = 0;
        condition = input_array[i];
        while(condition > 0){
          condition = input_array[i] -  pivots[j];
          j++;
        }
        bin = j-1;
        pos = dspls[bin] + scounts[bin] - bin_elements[bin];
        bucketlist[pos] = input_array[i];
        bin_elements[bin]--;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Alltoall(scounts,1,MPI_INT,local_size_array, 1, MPI_INT, MPI_COMM_WORLD);
//scatter array sizes
    int local_size = 0;
    for(i = 0 ; i < p ; i++){
      local_size += local_size_array[i];
    }
    local_dspls[0] = 0;
    for(i = 0 ; i< p-1 ;i++){
        local_dspls[i+1] = local_dspls[i] + local_size_array[i];
    }
    local_array = malloc(local_size*sizeof(double));
    MPI_Alltoallv(bucketlist, scounts, dspls, MPI_DOUBLE, local_array, local_size_array, local_dspls, MPI_DOUBLE, MPI_COMM_WORLD);
    tend1= MPI_Wtime();

//sort array
    tstart2 = MPI_Wtime();
    qsort_dbls(local_array, local_size);

//Amalgamate results
    MPI_Gather(&local_size,1,MPI_INT, scounts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    dspls[0] = 0;
    for(i = 0 ; i< p-1 ;i++){
        dspls[i+1] = dspls[i] + scounts[i];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(local_array, local_size, MPI_DOUBLE, orderedlist, scounts, dspls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    tend2 = MPI_Wtime();

    if(my_rank==0){
        printf("flag n p T2Gen T2Sort\n");
        printf("Data: %d %d %f %f\n",n,p,tend1-tstart1,tend2-tstart2);
    }
    int flag = 0;
    for(i = 0 ; i < n-1 ; i++){
        if(orderedlist[i]-orderedlist[i+1]>0){flag = 1;}
    }
    if(!flag){printf("success\n");}
    MPI_Finalize();

}

