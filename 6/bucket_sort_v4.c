#include <stdio.h>
#include <stdlib.h>
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

    MPI_Init(NULL,NULL);
    int p, my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    double *input_array;   // the input array
    double *bucketlist;    //this array will contain input elements in order of the processors
                            //e.g elements of process 0 will be stored first, then elements of process 1, and so on
    double *local_array;    //This array will contain the elements in each process
    double *orderedlist;
    double tstart1,tstart2, tend1, tend2;
    int *local_size_array;

    int n,i;

    int *scounts;           //This array will contain the counts of elements each processor will receive
    int *dspls;             //The relative offsets in bucketlist array where the elements of different processes
                            //will be stored
    int *local_dspls;
    int *bin_elements;      //it will keep track of how many elements have been included in the pth bi 
    n = atoi(argv[1]);

    tstart1 = MPI_Wtime();

//Memory Allocation
    int stop = n/p;
    if( my_rank == 0){stop += n%p;}
    input_array = malloc(stop*sizeof(double));
    orderedlist = malloc(n*sizeof(double));
    bucketlist = malloc(stop*sizeof(double));
    scounts = malloc(p*sizeof(int));
    local_size_array = malloc(p*sizeof(int));
    dspls = malloc(p*sizeof(int));
    local_dspls = malloc(p*sizeof(int));
    bin_elements = malloc(p*sizeof(int));

//Initialize each process random numbers
    double x=0;
    for(i = 0 ; i < stop ; i++){
        x = ((double) rand()/RAND_MAX);
        input_array[i] = x*x;
    }

    for(i = 0 ; i < p ; i++){
        scounts[i] = 0 ;
    }

//counting the elements in each processor
    for(i = 0 ; i < stop ; i++){
        scounts[(int)(input_array[i]*p)]++;
    }

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
        bin = (int)(input_array[i]*p);
        pos = dspls[bin] + scounts[bin] - bin_elements[bin];
        bucketlist[pos] = input_array[i];
        bin_elements[bin]--;
    }

//inform other processes how many elements are coming
    MPI_Alltoall(scounts,1,MPI_INT,local_size_array, 1, MPI_INT, MPI_COMM_WORLD);

//build array to hold incoming data
    int local_size = 0;
    for(i = 0 ; i < p ; i++){
      local_size += local_size_array[i];
    }
    local_dspls[0] = 0;
    for(i = 0 ; i< p-1 ;i++){
        local_dspls[i+1] = local_dspls[i] + local_size_array[i];
    }
    local_array = malloc(local_size*sizeof(double));

//Swap data to be sorted
    MPI_Alltoallv(bucketlist, scounts, dspls, MPI_DOUBLE, local_array, local_size_array, local_dspls, MPI_DOUBLE, MPI_COMM_WORLD);
    tend1 = MPI_Wtime();


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

//Check to make sure sorted correctly
    if(my_rank==0){
        printf("flag n p T2Gen T2Sort\n");
        printf("Data: %d %d %f %f\n",n,p,tend1-tstart1,tend2-tstart2);
        int flag = 0;
        for(i = 0 ; i < n-1 ; i++){
            if(bucketlist[i]-bucketlist[i+1]>0){flag = 1;}
        }
        if(!flag){printf("success\n");}
    }
    MPI_Finalize();

}

