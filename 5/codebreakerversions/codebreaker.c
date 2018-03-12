#include <inttypes.h> 
#include <stdlib.h> 
#include <stdio.h> 
#include <string.h> 
#include <math.h> 
#include <mpi.h>

#define ROTL32(x,y) ((x<<y)|(x>>(32-y)))
#define ROTR32(x,y) ((x>>y)|(x<<(32-y)))
#define ROTL24(x,y) ((x<<y)|(x>>(24-y)))
#define ROTR24(x,y) ((x>>y)|(x<<(24-y)))
#define ROTL16(x,y) ((x<<y)|(x>>(16-y)))
#define ROTR16(x,y) ((x>>y)|(x<<(16-y)))
#define ROTL8(x,y)  ((x<<y)|(x>>(8-y)))
#define ROTR8(x,y)  ((x>>y)|(x<<(8-y)))

int num_entries = 8;
const char * dictionary [] = {  "bird",
                                "campus",
                                "class",
                                "of",
                                "spring",
                                "sun",
                                "the",
                                "tree" };

// checks if a significant portion of the words resulting from
// the decryption with the tried key are found in the dictionary
int isValid( char *decoded, int len ) {

    int nmatches = 0;
    char * word;

    word = strtok( decoded, " ,.;-()\n\r" );

    while ( word != NULL ) {

        int flag = 0;

        for ( int i = 0; i < num_entries; ++i ) {
            if( !strcmp( word, dictionary[i] ) ) {
                flag = 1;
                break;
            }
        }

        if (flag) {  nmatches += strlen(word);  }

        word = strtok( NULL, " ,.;-()\n\r" );
    }

    // different criteria may be used for deciding whether the tried
    // key was a valid one. here we identify it as valid if words in
    // the decrypted message that can be located in the dictionary account
    // for more than half of the original message length.
    if (nmatches > len * 0.50)
        return 1;

    return 0;
}

void decrypt32( unsigned char *inp, uint32_t key, unsigned char *decoded ) {

    int i, iend, oend;
    uint32_t block;
    uint32_t a, b, c, d, magnitude, polarity, xor;

    srand(key);

    iend = 0;
    decoded[0] = 0; // C strings are zero-terminated
    oend = 0;

    /* main loop for decoding -- all 4 bytes are valid */
    while (     (a = inp[iend++]) != 0
            &&  (b = inp[iend++]) != 0
            &&  (c = inp[iend++]) != 0
            &&  (d = inp[iend++]) != 0
    ) {
        // printf("a = %x, b = %x, c = %x, d=%x\n", a, b, c, d);

        polarity = rand() % 2;
        magnitude = rand() % 32;

        block = (d << 24) | (c << 16) | (b << 8) | a;

        if (polarity)
            block = ROTR32(block, magnitude);
        else
            block = ROTL32(block, magnitude);

        xor =   (rand() % 256 << 24) | (rand() % 256 << 16)
              | (rand() % 256 << 8)  | (rand() % 256);

        block ^= xor;

        decoded[oend++] = block;
        decoded[oend++] = (block = block >> 8);
        decoded[oend++] = (block = block >> 8);
        decoded[oend++] = (block = block >> 8);
        decoded[oend] = 0;

        // printf("p = %d, mag = %d, xor = %d\n", polarity, magnitude, xor);
    }

    /* end cases */
    if ( a != 0 && b != 0 && c != 0 && d == 0 ) {

        polarity = rand() % 2;
        magnitude = rand() % 24;
        block = (c << 16) | (b << 8) | a;

        if (polarity)
            block = ROTR24(block, magnitude);
        else
            block = ROTL24(block, magnitude);

        xor = (rand() % 256 << 16) | (rand() % 256 << 8) | (rand() % 256);

        block ^= xor;

        decoded[oend++] = block;
        decoded[oend++] = (block = block >> 8);
        decoded[oend++] = (block = block >> 8);
        decoded[oend] = 0;

    } else if ( a != 0 && b != 0 && c == 0 ) {

        polarity = rand() % 2;
        magnitude = rand() % 16;
        block = (b << 8) | a;

        if (polarity)
            block = ROTR16(block, magnitude);
        else
            block = ROTL16(block, magnitude);

        xor = (rand() % 256 << 8) | (rand() % 256);

        block ^= xor;

        decoded[oend++] = block;
        decoded[oend++] = (block = block >> 8);
        decoded[oend] = 0;

    } else if ( a != 0 && b == 0 ) {

        polarity = rand() % 2;
        magnitude = rand() % 8;
        block = a;

        if (polarity)
            block = ROTR8(block, magnitude);
        else
            block = ROTL8(block, magnitude);

        xor = rand() % 256;
        block ^= xor;

        decoded[oend++] = block;
        decoded[oend] = 0;
    }
}


int main( int argc,char *argv[] ) {

    char a;
    char outfilename[100];
    unsigned char encrypted[1000], decrypted[1000], dcopy[1000];
    FILE *fin, *fout;
    int success = 0;
    uint64_t i, len, Nmax;
    double tstart, tend;

//    printf("\n\n x0r 32-bit code breaker\n\n");

    if (argc == 1) {
        fprintf(stderr, "ERROR: No file(s) supplied.\n");
        fprintf(stderr, "USAGE: This program requires a filename \
                        to be provided as argument!");
        exit(1);
    }

//    printf("decrypting file %s by trying all possible keys...\n", argv[1]);
//    printf("To quit, press ctrl + c\n\n");
//    printf("Status:\n");


    tstart = MPI_Wtime();

//Create Parallel environment
//MPI initialize Processes
   int Nprocess, Rprocess;
   MPI_Init(&argc, &argv);
    
   MPI_Comm_size(MPI_COMM_WORLD, &Nprocess);
   MPI_Comm_rank(MPI_COMM_WORLD, &Rprocess);
//   printf("From process %d out of %d, Hello World! \n", Rprocess, Nprocess);

//MPI initialize Communication style
   MPI_Request request;
   MPI_Status status;
   int flag=0;
   int index;

   MPI_Barrier(MPI_COMM_WORLD);

//process 0 reads in file and broadcasts
if(Rprocess==0){
   if ( (fin = fopen(argv[1], "r")) == NULL ) {
       fprintf(stderr, "ERROR: Could not open: %s\n", argv[1]);
       exit(1);
    }
    encrypted[0] = 0;
    len = 0;
    while ( (a = fgetc(fin)) != EOF )
        encrypted[len++] = a;
    encrypted[len] = 0;
    fclose(fin);
    MPI_Bcast(&encrypted, len, MPI_INT, Rprocess, MPI_COMM_WORLD);
}

MPI_Barrier(MPI_COMM_WORLD);

//each process launches pigeon post
printf("creating process: %d pigeon post \n", Rprocess);
MPI_Irecv(&i , 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);

//Define how frequently you want the processes to check pigeon post
#ifndef CHECK_INT
#define CHECK_INT 1000000
#endif

//    printf("starting decoder \n");
    Nmax = pow( 2, sizeof(uint32_t)*8);

    for (i=Nmax*Rprocess/Nprocess; i<Nmax*(Rprocess+1)/Nprocess; ++i) {

        decrypt32( encrypted, i, decrypted );

        strcpy(dcopy, decrypted);

        if (isValid(dcopy, len)) {
              success = 1;
              //send pigeons to all other posts to stop their process
              for(int k= 0;k<Nprocess;k++){if(k!=Rprocess)MPI_Send(&i, 1, MPI_INT, k, 123, MPI_COMM_WORLD);}
              break;
        }
        if((i+1)%CHECK_INT==0){
           printf("Rprocess : %d  is checking for mail \n",Rprocess);
           for(int k= 0;k<Nprocess;k++){if(k!=Rprocess){MPI_Test(&request, &flag, &status);}}
           if(flag){
              printf("Recieved carrier pigeon \n");
              MPI_Finalize();
              return 0;}
        }
    }

//END of MPI section 

    tend = MPI_Wtime();
    printf("\nTime elapsed: %.6e seconds\n", tend - tstart);

    if (success) {
        sprintf(outfilename, "%s.out", argv[1]);
        fout = fopen(outfilename, "w");
        fprintf(fout, "%s", decrypted);
        printf("\nFile decrypted successfully using key %u\n", (unsigned) i);
        printf("See the file %s\n\n\n", outfilename);
        fclose(fout);
        MPI_Finalize();
        return 0;
    }


    printf("\nWARNING: File could not be decrypted.\n\n\n");
    return 1;
}
