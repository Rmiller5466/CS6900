// cs6900
// Final - prob1_mpi
// Ryan Miller
// w051rem
// Due 27 April 2021
// System = bender
// Compiler syntax = ./prob1_mpi.compile prob1_mpi
// Job Control File = prob1_mpi.batch
// Additional File  = N/A
// Results file     = prob1_mpi.txt
// NOTE: Only the mersenne loop is parallelized.  

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <stdbool.h>
# include <time.h>
#include <mpi.h>

// ref http://www.tsm-resources.com/alists/mers.html

// Set below to 100 for testing your code
// When ready to run set this to 1 million 1000000
// submit one code called prob1_mpi.c by modifying
//        the code to use MPI and run faster
//        All cores must do work including rank zero
//        buffer tasks for each core
//        Use 2 nodes and all cores on the node.
// submit one called prob1_omp.c by modifying
//        the code to use openmp and run faster
//        Use 1 node and all cores on the node
// Note that your may only modify the 2 identified blocks
// All other parts of the code must be run by all ranks for MPI
// and a single thread for openMP.


#define MAXPRIME 100
#define MAXK MAXPRIME-2

int main ( int argc, char **argv );
bool is_prime( int n);
void make_prime_vector(int n, int *prime, int *k);
bool quick_is_prime(unsigned long long int j, int *prime, int k);

int main ( int argc, char **argv )
{
  int rank, ncpu;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ncpu);

  int prime[MAXPRIME];
  int mersenne[64];

  // Initial Prime values
  prime[0]=2;
  prime[1]=3;
  prime[2]=5;
  prime[3]=7;
  prime[4]=11;
  prime[5]=13;
  int k=6;
  int i;
  unsigned long long int j;

  clock_t begin = clock();

  // try to make the below to run in parallel
  // Create prime vector
  int n=17; // starting prime - skip even and factor of 6

  while (k<MAXK){
    make_prime_vector(n, prime, &k);
    make_prime_vector(n+2, prime, &k);
    n=n+6;
    /* The below many be helpful for debugging only
    if ((n-17)%10000==0)
      printf("n=%d k=%d\n",n,k);
    */
  }

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    // Only one core/thread prints the timing
  if (rank == 0){
    printf("time creating prime vector %f \n",time_spent);
  }

  // opmp has get see
  // https://www.openmp.org/spec-html/5.0/openmpsu160.html
  //try to modify to run in parallel
  n=0;

  int extra = 64 % ncpu;
  int split = 64 / ncpu;
  int st = 0;
  int en = 0;

  if (rank < extra) {
    st = 2 + rank * split ;
    en = st + split + 1;

  }else{
    st = 2 + rank * split;
    en = st + split;
  }
  
  int local[10] = {0};
  int localCount = 0;
  int globalCount = 0;
  int globalTemp[64] = {0};
  int currentCount = 0;
  int l;

  for (i=st; i<en; ++i){
    j=(unsigned long long int)pow(2,i)-1;
    if(quick_is_prime(j,prime,k)){
      local[localCount]=i;
      ++localCount;
    }
  }
  

  for (i = 0; i < ncpu; i++){
    if (i == rank){
      for (l = 0; l < localCount; l++){
	globalTemp[currentCount] = local[l];
	currentCount++;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&currentCount, 1, MPI_INT, i, MPI_COMM_WORLD);
  }

  MPI_Reduce(globalTemp, mersenne, 64, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&localCount, &globalCount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  n = globalCount;

  
    clock_t end2 = clock();
    double time_spent2 = (double)(end2 - end) / CLOCKS_PER_SEC;
    // Only one core/thread prints the timing
    if (rank == 0) {
      printf("time creating mersenne primes %f \n",time_spent2);
  }
  //output
  // Only one core/thread prints this output
  if (rank == 0){
    
    for (i=0; i<n; ++i){
      j=(unsigned long long int)pow(2,mersenne[i])-1;
      printf("2^(%d)-1 = %llu \n",mersenne[i],j);
    }
  }


  // Don't parallize this
  // It is not timed and very HARD
  int isum;
  long int sum=0;
  for (isum=0; isum<k; ++isum){
    if ( sum>1000000000)
      sum=sum-prime[isum];
    else
      sum=sum+prime[isum];

  }

  MPI_Barrier(MPI_COMM_WORLD);
  // All cores print this.... for openmp pnly one thread needs to print this
  printf("[%d] prime[%d]=%d sum=%d\n",rank,k-1,prime[k-1],sum);

  MPI_Finalize(); 
  return 0;
}

bool is_prime(int n){
  if (n <= 3)
    return(n > 1);
  if (n%2 == 0 || n%3 == 0)
    return(false);
  int i=5;
  while (i * i <= n){
    if (n%i == 0 || n%(i + 2) == 0)
      return(false);
    i=i+6;
  }
  return (true);
}

void make_prime_vector(int n, int *prime, int *k){

  int i=2;
  while ( i<(*k) ){
    if (n%prime[i]==0)
      return;
    ++i;
  }
  prime[(*k)]=n;
  ++(*k);
  return;
}


bool quick_is_prime(unsigned long long int j, int *prime, int k){

  int i=1;
  while ( i<k ){
    if (j%(unsigned long long int)prime[i]==0)
      return(false);
    if ((unsigned long long int)prime[i]*(unsigned long long int)prime[i]>j){
      return(true);
    }
    ++i;
  }

  unsigned long long int ii;
  ii=(unsigned long long int)(prime[k-2]+6);
  while (ii * ii <= j){
    if (j%ii == 0 || j%(ii + 2) == 0)
      return(false);
    ii=ii+6;
  }
  return (true);


}
