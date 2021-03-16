#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

// Flag 0: Random 1: Diagonal 2: Filename
// filename default is 0
void createMatrix(int N, double** mspace, double* ref, int flag, char* filename){ 

  int i, j;
  double sum; 

  if (flag == 0) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
	mspace[i][j] = (double)(rand() % 100);  
      }
    }
  }else if (flag == 1){
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
	if ( i == j) {
	  mspace[i][j] =  i + 1;
	}else{
	  mspace[i][j] = 0;
	}
      }
    }
  }else {
    mspace[0][0] = 2;
    mspace[0][1] = 3;
    mspace[0][2] = 1;
    mspace[0][3] = -1;
    mspace[1][0] = 6;
    mspace[1][1] = 0;
    mspace[1][2] = 15;
    mspace[1][3] = -6;
    mspace[2][0] = 4;
    mspace[2][1] = 3;
    mspace[2][2] = 5;
    mspace[2][3] = -4;
    mspace[3][0] = 8;
    mspace[3][1] = 18;
    mspace[3][2] = -2;
    mspace[3][3] = 3;
  }

  for (i = 0; i < N; i++){
    sum = 0;
    for (j = 0; j < N; j++){
      sum += mspace[i][j] * ref[j];
    }
    mspace[i][N] = sum;
  }
}

void printMat(double** mat, int N){
  int i,j;

  printf("Matrix %d\n", N);
  for (i = 0; i < N; i++) {
    for (j = 0; j < N+1; j++) {
      if (j == (N)){
	printf("%.3f\n", mat[i][j]);
      }else{
	printf("%.3f ", mat[i][j]);
      }
    }
  }
}

void printVec(double* vec, int N){
  int i;

  printf("Vector %d\n", N);
  for (i = 0; i < N; i++) {
      if (i == (N - 1)){
	printf("%.3f\n", vec[i]);
      }else{
	printf("%.3f ", vec[i]);
      }
  }
}

void GE(double** mat, int N){
  int base, row, col, swap;
  double maxVal;
  int maxIndex;
  double* temp;
  double cancelVal;

  for (base = 0; base < N-1; base++) {
    maxVal = fabs(mat[base][base]);
    maxIndex = base;
    for (swap = base; swap < N; swap++){
      if (fabs(mat[swap][base]) > maxVal){
    	maxVal = fabs(mat[swap][base]);
    	maxIndex = swap;
      }
    }
    if (maxIndex != base){
      temp = mat[base];
      mat[base] = mat[maxIndex];
      mat[maxIndex] = temp;
    }

    for (row = base + 1; row < N; row++){
      cancelVal = mat[row][base]/mat[base][base];
      for (col = 0; col < N+1; col++){
	mat[row][col] = mat[row][col] + (-1 * cancelVal  * mat[base][col]);
      }
    }
  }
}

void BS(double** mat, double* xVec, int N){
  int i, j;
  
  for ( i = N-1; i >= 0; i--){
    xVec[i] = mat[i][N] / mat[i][i];
    for (j = 0; j < i; j++){
      mat[j][N] = mat[j][N] - xVec[i] * mat[j][i];
      mat[j][i] = 0;
    }
  }
}



double main(double argc, char** argv) {
  srand(time(NULL));

  if (argc != 3) {
    fprintf(stderr, "Incorrect number of Inputs!\n");
    fprintf(stderr, "Usage: proj05S.exe -flag n\n");
    fprintf(stderr, "      -r n : Matrix A will use random numbers\n");
    fprintf(stderr, "      -d n : Matrix A will be a diagonal matrix\n");
    fprintf(stderr, "      -f filename : Read in Matrix A from a given file\n");
    //    exit(1);
  }

  char* nu = 'a';
  int N = 4;
  double* space1[N];
  double* ref1 = (double *)malloc(sizeof(double) * N);
  double* x1 = (double *)malloc(sizeof(double) * N);
  assert(ref1 != NULL);
  assert(x1 != NULL);
  int i, j;

  for (i = 0; i < N; i++){
    space1[i] = (double *)malloc(sizeof(double) * (N + 1));
    assert(space1[i] != NULL);
  }

  for (i = 0; i < N; i++) {
    ref1[i] = (double)(rand() % 15);  
  }

  /* ref1[0] = -27.0/2.0; */
  /* ref1[1] = 26.0/3.0; */
  /* ref1[2] = 7; */
  /* ref1[3] = 2; */

  printf("Reference Vector: \n");
  printVec(ref1, N);
  createMatrix(N, space1, ref1, 0, nu);
  printMat(space1, N);

  
  GE(space1, N);
  printf("After GE\n");
  printMat(space1, N);

  BS(space1, x1, N);
  printf("After BS\n");
  printf("X Vector:\n");
  printVec(x1, N);
  printMat(space1, N);


  for (i = 0; i < N; i++) {
    free(space1[i]);
    free(ref1);
    free(x1);
  }
}
