// cs6900
// Project 04a
// Ryan Miller
// Basis functions and matrix math functions based off of
//  examples Provided by: Dr. John Nehrbass
// w051rem
// Due 1 March 2021
// System = bender
// Compiler syntax = ./compileScalar.bash (or gcc proj04a.c -lm -o proj04a.exe)
// Job Control File = proj04a.sbatch
// Additional File  = compileScalar.bash
// Results file     = proj04a.txt


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>


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
//        l1 - int representing the number of rows
//        l2 - int representing the length of the rows
// Return: sums[2] - double vector
//                   sums[0]: Non-Diagonal Sum
//                   sums[1]: Diagonal Sum

double* calcMatrix(double **matrix, int l1, int l2){
 static double sums[2] = {0};
 int i,j;

 for (i = 0; i < l2; i++){
   for (j = 0; j < l1; j++){
     if (j == i){
       sums[1] += matrix[i][j];
     }else{
       sums[0] += matrix[i][j];
     }
      
   }
 }
return sums;
}

 


int main(int argc, char** argv) {


  if (argc != 7) {
    fprintf(stderr, "%d\n",argc);
    fprintf(stderr, "Usage: proj04a.exe L1 L2 A B Ap Bp\n");
    fprintf(stderr, "      L1,L2 :int - Length of row vectors 1 and 2\n");
    fprintf(stderr, "      A,B :int - Basis function coefficients for A and B\n");
    fprintf(stderr, "      Ap,Bp :int - Phase offset in degrees [0 360)\n");
    exit(1);
  }

  // Inputs
  int L1 = atoi(argv[1]);
  int L2 = atoi(argv[2]);
  int A  = atoi(argv[3]);
  int B  = atoi(argv[4]);
  int Ap = atoi(argv[5]);
  int Bp = atoi(argv[6]);


  // Constants
  double dx, dy, pi, norm1, norm2, a, b;
  pi = acos(-1.0);
  dx = (2.0*pi)/L1;
  dy = (2.0*pi)/L2;
  norm1 = sqrt(2.0/L1);
  norm2 = sqrt(2.0/L2);
  a = pi * Ap / 180.0;
  b = pi * Bp /180.0;
  
  // Print out constants for the run
  //printf("Pi: %f Dx: %f Dy: %f norm1: %f norm2: %f a: %f b: %f\n", pi, dx, dy, norm1, norm2, a, b);
  
  double *v1 = (double *)malloc(sizeof(double) * L1);
  assert(v1 != NULL);
  double *v2 = (double *)malloc(sizeof(double) * L1);
  assert(v2 != NULL);
  double *v3 = (double *)malloc(sizeof(double) * L2);
  assert(v3 != NULL);
  
  int i, j;
  double Xi, Yi;

  for (i = 0; i < L1; i++) {
    Xi = dx/2.0 + i * dx;
    v1[i] = norm1*cos((A*Xi)+a);
    v2[i] = norm1*sin((B*Xi)+b);
  }

  for (i = 0; i < L2; i++){
    Yi = dy/2.0 + i * dy;
    v3[i] = norm2*cos((B*Yi)+b);
  }

  // Print out the vectors
  
  /*  printf("V1: "); */
  /*  for (i = 0; i < L1; i++) { */
  /*    printf("%f ",v1[i]); */
  /*  } */
  /*  printf("\n"); */
  
  /* printf("V2: "); */
  /*  for (i = 0; i < L1; i++) { */
  /*    printf("%f ",v2[i]); */
  /*  } */
  /*  printf("\n"); */
  
  /* printf("V3: "); */
  /*  for (i = 0; i < L2; i++) { */
  /*    printf("%f ",v3[i]); */
  /*  } */
  /*  printf("\n"); */
  
  double dotSum;
  dotSum = dot(v1,v2,L1);

  double *Matrix[L2];
  assert(Matrix != NULL);

  for (i = 0; i < L2; i++){
    Matrix[i] = (double *)malloc(sizeof(double) * L1);
    assert(Matrix[i] != NULL);
  }

  for (i = 0; i < L2; i++){
    for (j = 0; j < L1; j++){
      Matrix[i][j] = v3[i] * v1[j];
    }
  }

  double *matrixTotals = calcMatrix(Matrix, L1, L2);

  printf("L1: %d L2: %d A: %d B: %d Ap: %d Bp: %d\n",L1, L2, A, B, Ap, Bp);
  printf("The dot product is %.30f\n",dotSum);
  printf("Matrix Total: %.30f Matrix Diagonal Total: %.30f\n\n",matrixTotals[0], matrixTotals[1]);
  
  // Clean up Memory
  free(v1);
  free(v2);
  free(v3);
  for (i = 0; i < L2; i++){
    free(Matrix[i]);
  }
}
