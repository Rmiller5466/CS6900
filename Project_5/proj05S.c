#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

// Flag 0: Random 1: Diagonal 2: Filename
// filename default is 0
void createMatrix(int N, double** mspace, int flag, char* filename){ 

  int i, j;

  if (flag == 0) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
	mspace[i][j] = (double)(rand() % 50);  
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
    mspace[0][0] = 1;
    mspace[0][1] = 2;
    mspace[0][2] = 3;
    mspace[1][0] = 2;
    mspace[1][1] = -1;
    mspace[1][2] = 1;
    mspace[2][0] = 3;
    mspace[2][1] = 0;
    mspace[2][2] = -1;





  }
}

void printMat(double** mat, int N){
  int i,j;

  printf("Matrix %d\n", N);
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      if (j == (N - 1)){
	printf("%.2f\n", mat[i][j]);
      }else{
	printf("%.2f ", mat[i][j]);
      }
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
    /* maxVal = mat[base][base]; */
    /* maxIndex = base; */
    /* for (swap = base; swap < N; swap++){ */
    /*   if (mat[swap][base] > maxVal){ */
    /* 	maxVal = mat[swap][base]; */
    /* 	maxIndex = swap; */
    /*   } */
    /* } */
    /* if (maxIndex != base){ */
    /*   temp = mat[base]; */
    /*   mat[base] = mat[maxIndex]; */
    /*   mat[maxIndex] = temp; */
    /* } */

    /* for (j = i + 1; j < N-1; j++) { */
    /*   temp = mat[j][i]/mat[i][i]; */
    /*   for (k = i+1; j < N-1; j++) { */
    /* 	printf("I: %d, J: %d, K: %d\n",i,j,k); */
    /* 	mat[j][k] = ( mat[j][k] - (temp * mat[i][k])); */
    /*   }  */
    /* } */
    for (row = base + 1; row < N; row++){
      for (col = 0; col < N; col++){
	mat[row][col] = mat[row][col] + (-1 * (mat[row][base]/mat[base][base]) * mat[base][col]);
      }
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
  int N = 6;
  double* space1[N];
  int i, j;

  for (i = 0; i < N; i++){
    space1[i] = (double *)malloc(sizeof(double) * N);
    assert(space1[i] != NULL);
  }

 
  createMatrix(N, space1, 0, nu);
  printMat(space1, N);
  GE(space1, N);
  printMat(space1, N);
  for (i = 0; i < N; i++) {
    free(space1[i]);
  }

}
