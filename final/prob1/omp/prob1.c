# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <stdbool.h>
# include <time.h>
# include <omp.h>

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
  int prime[MAXPRIME];
  int mersenne[64];
  int rank = 0;

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
  while ( k<MAXK){
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
  printf("time creating prime vector %f \n",time_spent);

  //  clock_t end = clock();
  // opmp has get see
  // https://www.openmp.org/spec-html/5.0/openmpsu160.html
  //try to modify to run in parallel
  n=0;
  for (i=2; i<64; ++i){
    j=(unsigned long long int)pow(2,i)-1;
    if(quick_is_prime(j,prime,k)){
      mersenne[n]=i;
      ++n;
    }
  }

  clock_t end2 = clock();
  double time_spent2 = (double)(end2 - end) / CLOCKS_PER_SEC;
  // Only one core/thread prints the timing
  printf("time creating mersenne primes %f \n",time_spent2);

  //output
  // Only one core/thread prints this output
  for (i=0; i<n; ++i){
    j=(unsigned long long int)pow(2,mersenne[i])-1;
    printf("2^(%d)-1 = %llu \n",mersenne[i],j);
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

  // All cores print this.... for openmp pnly one thread needs to print this
  printf("[%d] prime[%d]=%d sum=%d\n",rank,k-1,prime[k-1],sum);

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
