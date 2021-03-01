// cs6900
// Project 04b
// Ryan Miller
// Basis functions and matrix math functions based off of
//  examples Provided by: Dr. John Nehrbass
// w051rem
// Due 1 March 2021
// System = bender
// Compiler syntax = ./compile.sh proj04b
// Job Control File = proj04b.sbatch
// Additional File  = NA
// Results file     = proj04b.txt

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

// dot: Calculates the dot product of two vectors
// Input: v1 - vector of doubles
//        v2 - vector of doubles
//        N - length of both vectors
// Return: sum - double representing the dot product

double dot(double *v1, double *v2, int N) {
  int i;
  double sum=0.0;
  for (i = 0; i < N; i++) {
    double temp = v1[i]*v2[i];
    sum += v1[i]*v2[i];
  }
  return sum;
}

// calcMatrix: Calculates the matrix formed by multiplying
//             two vectors
// Input: matrix - mxn vector of doubles 
//        Nr - int representing the number of rows
//        l2 - int representing the length of the rows
//        offset - used to calculate the diagonal 
//        sums - pointer which will contain the results 

void calcMatrix(double **matrix, int Nr, int l2, int offset, double *sums){
  int i,j;

  for (i = 0; i < Nr; i++){
    for (j = 0; j < l2; j++){
      if (j == (offset-1+i)){
	sums[1] += matrix[i][j];
      }else{
	sums[0] += matrix[i][j];
      }
    }
  }
}

 


int main(int argc, char** argv) {

  int rank, ncpu;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ncpu);
  
  // The next two blocks can be used to read in input. 
  // Note: Further manipulation would be needed to make it work
  //       (Taking out the following for loop)

  /* if (argc != 7 && rank == 0) { */
  /*   fprintf(stderr, "%d\n",argc); */
  /*   fprintf(stderr, "Usage: proj04a.exe L1 L2 A B Ap Bp\n"); */
  /*   fprintf(stderr, "      L1,L2 :int - Length of row vectors 1 and 2\n"); */
  /*   fprintf(stderr, "      A,B :int - Basis function coefficients for A and B\n"); */
  /*   fprintf(stderr, "      Ap,Bp :int - Phase offset in degrees [0 360)\n"); */
  /*   MPI_Abort(MPI_COMM_WORLD,1); */
  /*   exit(1); */
  /* } */

  /* // Inputs */
  /* int L1 = atoi(argv[1]); */
  /* int L2 = atoi(argv[2]); */
  /* int A  = atoi(argv[3]); */
  /* int B  = atoi(argv[4]); */
  /* int Ap = atoi(argv[5]); */
  /* int Bp = atoi(argv[6]); */

  // 10^1 - 10^9
  int L1_CONST[] =
    {10,
     100,
     1000,
     10000,
     100000,
     1000000,
     10000000,
     100000000,
     1000000000};
  
  int MAT_CONST[] = {1000, 100000};
  int OFFSET_CONST[] = {22, 45, 75, 90, 135};
  int runs;
  
  int L1, L2, A, Ap, Bp;
  long B = 0;
  L2 = 0;

  // Controls inputs for each iteration of the program

  for (runs = 0; runs < 31; runs++){
    if (runs < 9){

      L1 = L1_CONST[runs % 9];
      A = 1; B = -1; Ap = OFFSET_CONST[1]; Bp = OFFSET_CONST[1];
    
    }
    if (runs > 8 && runs <18){
    
      L1 = L1_CONST[runs % 9];
      A = 1; B = (long)(L1/4); Ap = OFFSET_CONST[0]; Bp = OFFSET_CONST[3];

    }
    if (runs > 17 && runs < 27){
      
      L1 = L1_CONST[runs % 9];
      A = 1; B = (long)(2.2*L1); Ap = OFFSET_CONST[2]; Bp = OFFSET_CONST[4];
    
    }
    if (runs > 26 && runs < 29){
    
      L1 = MAT_CONST[runs % 2]; L2 = L1;
      A = 1; B = 1; Ap = OFFSET_CONST[4]; Bp = OFFSET_CONST[4];
    
    }
    if (runs > 28 && runs < 31){
      L1 = MAT_CONST[runs % 2]; L2 = L1;
      A = 1; B = (long)(2.4*L2); Ap = OFFSET_CONST[1]; Bp = OFFSET_CONST[4];
   
    }
    

  
    
    // Constants
    double dx, dy, pi, norm1, norm2, a, b;
    pi = acos(-1.0);
    dx = (2.0*pi)/L1;
    dy = (2.0*pi)/L2;
    norm1 = sqrt(2.0/L1);
    norm2 = sqrt(2.0/L2);
    a = pi * Ap / 180.0;
    b = pi * Bp /180.0;
  
  
    // Calculates which processor has what part of the global elements
    // Block row Partitioned
    int Nr, R, row;
    Nr= (int)(L1/ncpu);   // All processors have at least this many rows
    R = (L1%ncpu);       // Remainder of above calculation
    if (rank<R) {
      ++Nr;             // First R processors have one more element
      row=1+Nr*rank;    // starting row
    }else{
      row=1+(Nr+1)*R+Nr*(rank-R); // starting row
    }

    double *v1 = (double *)malloc(sizeof(double) * Nr);
    assert(v1 != NULL);
    double *v2 = (double *)malloc(sizeof(double) * Nr);
    assert(v2 != NULL);
    double *v3 = (double *)malloc(sizeof(double) * L2);
    assert(v3 != NULL);

    int i, j;
    double Xi, Yi;

    for (i = 0; i < Nr; i++)  {
      Xi = dx/2.0 + (row-1+i)  * dx;
      v1[i] = norm1*cos((A*Xi)+a);
      v2[i] = norm1*sin((B*Xi)+b);
    }

    for (i = 0; i < L2; i++){
      Yi = dy/2.0 + i * dy;
      v3[i] = norm2*cos((B*Yi)+b);
    }
    
    // Print out the vectors on each rank
  
    /* printf("Rank: %d, V1: ", rank); */
    /* for (i = 0; i < Nr; i++) { */
    /*   printf("%f ",v1[i]); */
    /* } */
    /* printf("\n"); */

    /* printf("Rank: %d, V2: ", rank); */
    /* for (i = 0; i < Nr; i++) { */
    /*   printf("%f ",v2[i]); */
    /* } */
    /* printf("\n"); */

    /* printf("Rank: %d, V3: ", rank); */
    /* for (i = 0; i < L2; i++) { */
    /*   printf("%f ",v3[i]); */
    /* } */
    /* printf("\n"); */
  
    double localDotSum;
    localDotSum = dot(v1,v2,Nr);
  
    double *Matrix[Nr];
    assert(Matrix != NULL);

    for (i = 0; i < Nr; i++){
      Matrix[i] = (double *)malloc(sizeof(double) * L2);
      assert(Matrix[i] != NULL);
    }

    for (i = 0; i < Nr; i++){
      for (j = 0; j < L2; j++){
	Matrix[i][j] = v3[j] * v1[i];
      }
    }
  
    // Can be used to see the matrix being generated
  
    /*   MPI_Barrier(MPI_COMM_WORLD); */
   
    /*   if( rank == 0) { */
    /*     printf("Rank: %d, Matrix: \n", rank); */
    /*     for (i = 0; i < Nr; i++) { */
    /*       for ( j = 0; j < L2; j++){ */
    /* 	 printf("%f ", Matrix[i][j]); */
    /*       } */
    /*       printf("\n"); */
    /*     } */
    /*     printf("\n"); */
    /*   } */
   
    /*   MPI_Barrier(MPI_COMM_WORLD); */


    double localMatrixTotals[2] = {0};
    calcMatrix(Matrix, Nr, L2, row, localMatrixTotals);

    double dotSumRecv;
    double matrixTotalsRecv[2] = {0};

    // Condense all totals down to rank 0 and print out results

    if (rank == 0){
      MPI_Reduce(&localDotSum,&dotSumRecv,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(localMatrixTotals,matrixTotalsRecv,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

      // Output constants used: 
      // printf("Pi: %f Dx: %f Dy: %f norm1: %f norm2: %f a: %f b: %f\n", pi, dx, dy, norm1, norm2, a, b);

      printf("L1: %d L2: %d A: %d B: %ld Ap: %d Bp: %d\n",L1, L2, A, B, Ap, Bp);
      printf("Dot Product of V1:V2 was: %.30f\n",dotSumRecv);
      printf("Matrix Total (Non-Diagonal): %.30f, Matrix Diagonal Total: %.30f\n\n",matrixTotalsRecv[0], matrixTotalsRecv[1]);
   
      MPI_Barrier(MPI_COMM_WORLD);
  
    }else{
      MPI_Reduce(&localDotSum,&dotSumRecv,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(localMatrixTotals,matrixTotalsRecv,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
    }
  
    MPI_Barrier(MPI_COMM_WORLD);
 
    //Clean up Memory
    free(v1);
    free(v2);
    free(v3);
    for (i = 0; i < Nr; i++){
      free(Matrix[i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  // For loop ended

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}
