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
    double tstart1,tstart2, tend1, tend2;
    int local_size;

    int n,i;

    int *scounts;           //This array will contain the counts of elements each processor will receive
    int *dspls;             //The relative offsets in bucketlist array where the elements of different processes
                            //will be stored
    int *bin_elements;      //it will keep track of how many elements have been included in the pth bi 
    n = atoi(argv[1]);

//Memory allocation
    tstart1 = MPI_Wtime();

    input_array = malloc(n*sizeof(double));
    bucketlist = malloc(n*sizeof(double));
    scounts = malloc(p*sizeof(int));
    dspls = malloc(p*sizeof(int));
    bin_elements = malloc(p*sizeof(int));

//Rank 0 gnerates the random numbers and distributes them
if(my_rank==0){
    //generate random data
    for(i = 0 ; i < n ; i++){
        input_array[i] = ((double) rand()/RAND_MAX);
    }

    //counting the elements in each processor
    for(i = 0 ; i < p ; i++){
        scounts[i] = 0 ;
    }
    for(i = 0 ; i < n ; i++){
        scounts[(int)(input_array[i]/(1.0/p))]++;
    }

    //Place elements into buckets
    for(i = 0 ; i<p ; i++){
        bin_elements[i] = scounts[i];
    }

    dspls[0] = 0;
    for(i = 0 ; i< p-1 ;i++){
        dspls[i+1] = dspls[i] + scounts[i];
    }
    int bin;
    int pos;
    for(i = 0 ; i < n ; i++){
        bin = (int)(input_array[i]/(1.0/p));
        pos = dspls[bin] + scounts[bin] - bin_elements[bin];
        bucketlist[pos] = input_array[i];
        bin_elements[bin]--;
    }
}


//scatter array sizes
    MPI_Scatter(scounts,1,MPI_INT,&local_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

//make array large enough to hold incoming data
    local_array = malloc(local_size*sizeof(double));

//rank 0 sends each process it bucket to sort
    MPI_Scatterv(bucketlist, scounts, dspls, MPI_DOUBLE, local_array, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    tend1=MPI_Wtime();

//sort array
    tstart2 = MPI_Wtime();
    qsort_dbls(local_array, local_size);

//Amalgamate results
    MPI_Gatherv(local_array, local_size, MPI_DOUBLE, bucketlist, scounts, dspls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    tend2 = MPI_Wtime();

//Check to make sure sort was done correctly
   if(my_rank==0){
   printf("flag n p T2Gen T2Sort TotalT\n");
   printf("Data: %d %d %f %f\n",n,p,tend1-tstart1,tend2-tstart2);
   int flag = 0;
   for(i = 0 ; i < n-1 ; i++){
       if(bucketlist[i]-bucketlist[i+1]>0){
           flag = 1;
       }

   }
   if(!flag){printf("success\n");}
   }

//close communication interface
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

}

