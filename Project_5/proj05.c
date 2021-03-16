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

void createMatrix(int N, int N_global, double** mspace, char* flag, char* filename, int rank, int ncpu){ 

  int i, j; // k;
  double sum; 
  int offset = rank;

  if (strcmp(flag, "-r") == 0) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < N_global; j++) {
	mspace[i][j] = (double)(rand() % 100);  
      }
    }
  }else if (strcmp(flag, "-d") == 0){
    for (i = 0; i < N; i++) {
      offset = rank + (i * ncpu);
      for (j = 0; j < N_global; j++) {
	if ( offset == j) {
	  mspace[i][j] =  offset + 1;
	}else{
	  mspace[i][j] = 3.141;
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
    
    int counter = 0;
      
    num = strtok(st, " ");
    for (i = 0; i < N_global; i++) { 
      for (j = 0; j < N_global; j++) {
	if (i == offset){
	  mspace[counter][j] = atof(num);	
	}
	if (i == offset && j == N_global - 1){
	  counter++;
	  offset = rank + (counter * ncpu);
	}
	num = strtok(NULL, " ");
      }
    }


  }else{
    printf("Emergency Abort: Invalid Input Flag (Must use -r, -d, or -f!)\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
  }
}

void createRefMatrix(int N, int M,  double** ref, char* flag, char* filename, int rank, int ncpu){ 
  int i, j, lines;
  int offset = rank;
	

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
    int counter = 0;

    fgets(st, sizeof(st), fp);
    for (i = 0; i < M; i++){
      eof = fgets(st,sizeof(st), fp);
      //      printf("Ref String: %s\n",st);
      
      if (eof == NULL) {
	printf("Emergency Abort: Input File Line Count Didn't Match M!\n");
	MPI_Abort(MPI_COMM_WORLD, 1);
	exit(1);
      }
      num = strtok(st, " ");
      
      for (j = 0; j < N; j++) {
	//	printf("Num: %s Offset:  Counter:  \n", num);// offset, counter);
	if (j == offset){
	  ref[counter][i] = atof(num);
	  counter++;
	  offset = rank + (counter * ncpu);
	}
	//if (j == offset &&
	//ref[j][i] = atof(num);  
	//ref[j][i] = atof(num);
	num = strtok(NULL, " ");
      }
      counter = 0;
      offset = rank;
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

  int rank, ncpu;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ncpu);

  srand(time(NULL)*(double)rank);
  
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

  // Calculates which processor has what part of the global elements
  int local_n, R;
  local_n= (int)(n/ncpu);   // All processors have at least this many rows
  R = (n%ncpu);       // Remainder of above calculation
  if (rank<R) {
    ++local_n;             // First R processors have one more element
  }

  double** A = (double**)malloc(local_n * sizeof(double*));
  double** r = (double**)malloc(local_n * sizeof(double*));
  double** x = (double**)malloc(local_n * sizeof(double*));
  int i, j;

  for (i = 0; i < local_n; i++){
    A[i] = (double *)malloc(sizeof(double) * (n + m));
    r[i] = (double *)malloc(sizeof(double) * m);
    x[i] = (double *)malloc(sizeof(double) * m);
    assert(A[i] != NULL);
    assert(r[i] != NULL);
    assert(x[i] != NULL);
  }

  /* for (i = 0; i < local_n; i++) { */
  /*   for (j = 0; j < m; j++){ */
  /*     r[i][j] = (double) (rand()%15); */
  /*   } */
  /* } */

  createMatrix(local_n, n, A, matFlag, fileNameMat, rank, ncpu);
  createRefMatrix(n, m, r, refFlag, fileNameRef, rank, ncpu);

  MPI_Barrier(MPI_COMM_WORLD);
   
  if( rank == 0) {
    
  printf("Rank %d's matrix:\n", rank);
  printMat(A, local_n, n+m);
  printMat(r, local_n, m);
  }

  MPI_Barrier(MPI_COMM_WORLD);
   
  if( rank == 1) {
    
  printf("Rank %d's matrix:\n", rank);
  printMat(A, local_n, n+m);
  printMat(r, local_n, m);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if( rank == 2) {
    
  printf("Rank %d's matrix:\n", rank);
  printMat(A, local_n, n+m);
  printMat(r, local_n, m);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  //printf("Rank %d's matrix:\n", rank);
  //printMat(A, local_n, n+m);
  /* printf("Rank %d's reference:\n", rank); */
  
  /* generateY(n, m, A, r); */
  /* printf("Rank %d's Full Matrix:\n", rank); */
  /* printMat(A,n, n+m);; */


  for (i = 0; i < local_n; i++) {
    free(A[i]);
    free(r[i]);
    free(x[i]);
  }
  free(A);
  free(r);
  free(x);

  MPI_Finalize(); 
}
