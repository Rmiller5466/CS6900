
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

float dot(float *v1, float *v2, int N) {
  int i;
  float sum=0.0;
  for (i = 0; i < N; i++) {
    sum += v1[i]*v2[i];
  }
  return sum;
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
  float dx, dy, pi, norm1, norm2, a, b;
  pi = acos(-1.0);
  dx = (2.0*pi)/L1;
  dy = (2.0*pi)/L2;
  norm1 = sqrt(2.0/L1);
  norm2 = sqrt(2.0/L2);
  a = pi*(Ap/180.0);
  b = pi*(Bp/180.0);

  float *v1 = (float *)malloc(sizeof(float) * L1);
  assert(v1 != NULL);
  float *v2 = (float *)malloc(sizeof(float) * L1);
  assert(v2 != NULL);
  float *v3 = (float *)malloc(sizeof(float) * L2);
  assert(v3 != NULL);
  
  int i;
  float Xi, Yi;
  
  for (i = 0; i < L1; i++) {
    Xi = Xi+dx;
    Yi = Yi+dy;
    
    v1[i] = norm1*cos((A*Xi)+a);
    v2[i] = norm1*sin((B*Xi)+b);
    // v3[i] = norm2*cos((B*Yi)+b);
  }

  // For loop L2

  float dotSum;
  dotSum = dot(v1,v2,L1);

  printf("L1: %d\n",L1);
  printf("The dot product is %f\n",dotSum);


  // Clean up Memory
  free(v1);
  free(v2);

}
