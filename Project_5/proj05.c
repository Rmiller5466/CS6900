// cs6900
// Project 05
// Ryan Miller
// w051rem
// Due 19 March 2021
// System = bender
// Compiler syntax = ./compile.sh proj05
// Job Control File = proj05.sbatch
// Additional File  = NA
// Results file     = proj05.txt

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <string.h>

void printVec(double* vec, int N){
  int i;

  for (i = 0; i < N; i++) {
    if (i == (N - 1)){
      printf("%.3f\n", vec[i]);
    }else{
      printf("%.3f ", vec[i]);
    }
  }
}

void printMat(double** mat, int N, int M){
  int i,j;

  printf("Matrix Length: %d\n", N);
  for (i = 0; i < N; i++) {
    for (j = 0; j < M; j++) {
      if (j == (M - 1)){
	printf("%.3f\n", mat[i][j]);
      }else{
	printf("%.3f ", mat[i][j]);
      }
    }
  }
}

// Flag 0: Random 1: Diagonal 2: Filename
// filename default is n/m
//void createMatrix(int N, int M, double** mspace, double** ref, int flag, char* filename){ 
void createMatrix(int N, double** mspace, char* flag, char* filename){ 

  int i, j; // k;
  double sum; 


  if (strcmp(flag, "-r") == 0) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
	mspace[i][j] = (double)(rand() % 100);  
      }
    }
  }else if (strcmp(flag, "-d") == 0){
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
	if ( i == j) {
	  mspace[i][j] =  i + 1;
	}else{
	  mspace[i][j] = 0;
	}
      }
    }
  }else if(strcmp(flag, "-f") == 0) {
    
    FILE* fp;
    
    fp = fopen(filename, "r"); 

    if (fp == NULL){
      printf("Emergency Abort: File Not Found!\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(1);
    }

    char st[4096];
    char* num;
    char* eof;

    fgets(st, sizeof(st), fp);
    eof = fgets(st,sizeof(st), fp);
    if (eof == NULL) {
      printf("Emergency Abort: Incorrect File Format!\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(1);
    }
    
    fclose(fp);
    printf("String: %s\n",st);
    
    num = strtok(st, " ");
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
	mspace[i][j] = atof(num);  
	num = strtok(NULL, " ");
      }
    }


  }else{
    printf("Emergency Abort: Invalid Input Flag (Must use -r, -d, or -f!)\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
  }
}

void createRefMatrix(int N, int M,  double** ref, char* flag, char* filename){ 
  int i, j, lines;

  if (strcmp(flag, "-R") == 0) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < M; j++) {
	ref[i][j] = (double)(rand() % 100);  
      }
    }
  }else if(strcmp(flag, "-F") == 0) {
    
    FILE* fp;
    
    fp = fopen(filename, "r"); 

    if (fp == NULL){
      printf("Emergency Abort: File Not Found!\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(1);
    }

    char st[4096];
    char* eof;
    char* num;
    fgets(st, sizeof(st), fp);
    for (i = 0; i < M; i++){
      eof = fgets(st,sizeof(st), fp);
      if (eof == NULL) {
	printf("Emergency Abort: Input File Line Count Didn't Match M!\n");
	MPI_Abort(MPI_COMM_WORLD, 1);
	exit(1);
      }
      num = strtok(st, " ");
      
      for (j = 0; j < N; j++) {
	ref[j][i] = atof(num);  
	num = strtok(NULL, " ");
      }
    }
    fclose(fp);

  }else{
    printf("Emergency Abort: Invalid Input Flag (Must use -R, or -F!)\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
  }
  
}

void generateY(int N, int M, double** mspace, double** ref){
  int i, j, k;
  double sum;
  
  for (k = 0; k < M; k++){
    for (i = 0; i < N; i++){
      sum = 0;
      for (j = 0; j < N; j++){
  	sum += mspace[i][j] * ref[i][k];
      }
      mspace[i][N+k] = sum;
    }
  }
}

int main(int argc, char** argv) {

  srand(time(NULL));

  int rank, ncpu;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ncpu);

  if (rank == 0 && argc != 5) {
    fprintf(stderr, "Incorrect number of Inputs!\n");
    fprintf(stderr, "Usage: proj05S.exe -MatrixFlag n -ReferenceFlag m\n");
    fprintf(stderr, "   -r n : Matrix A will use random numbers\n");
    fprintf(stderr, "   -d n : Matrix A will be a diagonal matrix\n");
    fprintf(stderr, "   -f filename : Read in Matrix A from a given file\n");
    fprintf(stderr, "   -R m : Reference Matrix of nxm random numbers\n");
    fprintf(stderr, "   -F filename : Read in Reference Matrix from a given file\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  char* matFlag;
  char* refFlag;
  char* fileNameMat;
  char* fileNameRef;
  int n = 0, m = 0;
  
  matFlag = argv[1];
  fileNameMat = argv[2];

  if (strcmp(matFlag, "-f") == 0){
    FILE* fp;
    char tmp[64];

    fp = fopen(fileNameMat, "r");
    if (fp == NULL){
      printf("Emergency Abort: File Not Found!\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(1);
    } 
    fgets(tmp, sizeof(tmp), fp);
    n = atoi(tmp);
    fclose(fp);
  }else{
    n = atoi(argv[2]);
  }
  
  refFlag = argv[3];
  fileNameRef = argv[4];
    
  if (strcmp(refFlag, "-F") == 0){
     FILE* fp;
     char tmp[64];
    fp = fopen(fileNameRef, "r");
    if (fp == NULL){
      printf("Emergency Abort: File Not Found!\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(1);
    } 
    fgets(tmp, sizeof(tmp), fp);
    m = atoi(tmp);
    
    fclose(fp);

  }else{
    m = atoi(argv[4]);
  }

  printf("Inputs:  Flag1: %s n: %d fileMat: %s Flag2: %s m: %d fileRef: %s\n", matFlag, n, fileNameMat, refFlag, m, fileNameRef);


  double** A = (double**)malloc(n * sizeof(double*));
  double** r = (double**)malloc(n * sizeof(double*));
  double** x = (double**)malloc(n * sizeof(double*));
  int i, j;

  for (i = 0; i < n; i++){
    A[i] = (double *)malloc(sizeof(double) * (n + m));
    r[i] = (double *)malloc(sizeof(double) * m);
    x[i] = (double *)malloc(sizeof(double) * m);
    assert(A[i] != NULL);
    assert(r[i] != NULL);
    assert(x[i] != NULL);
  }

  for (i = 0; i < n; i++) {
    for (j = 0; j < m; j++){
      r[i][j] = (double) (rand()%15);  
    }
  }

  createMatrix(n, A, matFlag, fileNameMat);
  createRefMatrix(n, m, r, refFlag, fileNameRef);
  printMat(A, n, n+m);
  printMat(r, n, m);
  generateY(n, m, A, r);
  printMat(A,n, n+m);;


  for (i = 0; i < n; i++) {
    free(A[i]);
    free(r[i]);
    free(x[i]);
  }
  free(A);
  free(r);
  free(x);

  MPI_Finalize(); 
}
