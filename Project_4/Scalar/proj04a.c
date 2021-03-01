#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

double dot(double *v1, double *v2, int N) {
  int i;
  double sum=0.0;
  for (i = 0; i < N; i++) {
    double temp = v1[i]*v2[i];
    //printf("i value: %d, %f ",i,temp);
    sum += v1[i]*v2[i];
  }
  // printf("\n\n");
  return sum;
}

double* calcMatrix(double **matrix, int l1, int l2){
 static double sums[2] = {0};
 int i,j;

 for (i = 0; i < l2; i++){
   for (j = 0; j < l1; j++){
     printf("%f ", matrix[i][j]);
     if (j == i){
       //printf("%f\n", matrix[i][j]);
       sums[1] += matrix[i][j];
     }else{
       sums[0] += matrix[i][j];
     }
      
   }
   printf("\n");
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
  
  printf("Pi: %f Dx: %f Dy: %f norm1: %f norm2: %f a: %f b: %f\n", pi, dx, dy, norm1, norm2, a, b);
  
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

  printf("V1: ");
  for (i = 0; i < L1; i++) {
    printf("%f ",v1[i]);
  }
  printf("\n");

 printf("V2: ");
  for (i = 0; i < L1; i++) {
    printf("%f ",v2[i]);
  }
  printf("\n");

 printf("V3: ");
  for (i = 0; i < L2; i++) {
    printf("%f ",v3[i]);
  }
  printf("\n");

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
      // printf("v1 J-%d VALUE: %f\nv3 I-%d VALUE: %f\nMATRIX-I: %d J: %d VALUE:%f\n\n",j,v1[j],i,v3[i],i,j,Matrix[i][j]);
    }
  }

  double *matrixTotals = calcMatrix(Matrix, L1, L2);

  printf("L1: %d L2:%d\n",L1,L2);
  printf("The dot product is %f\n",dotSum);
  printf("Matrix Total: %f Matrix Diagonal Total: %f\n\n",matrixTotals[0], matrixTotals[1]);
  
  // Clean up Memory
  free(v1);
  free(v2);
  free(v3);
  for (i = 0; i < L2; i++){
    free(Matrix[i]);
  }
}
