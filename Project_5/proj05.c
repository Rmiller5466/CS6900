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

#define DEBUG_PRINT_MAT
#define DEBUG_PRINT_RECV
#define DEBUG_INPUTS
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

void printMat(double** mat, int local_N, int N, int M){
  int i,j;

  for (i = 0; i < local_N; i++) {
    for (j = 0; j < N+M; j++) {
      if (j == (N+M - 1)){
  	printf("%.3f\n", mat[i][j]);
      }else{
  	if (j == N){
  	  printf(": ");
  	}
  	printf("%.3f ", mat[i][j]);
      }
    }
  }
}

void createMatrix(int N, int N_global, double** mspace, char* flag, char* filename, int rank, int ncpu){ 

  int i, j; 
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

void createRefMatrix(int l_n, int N, int M,  double** ref, char* flag, char* filename, int rank, int ncpu){ 
  int i, j, lines;
  int offset = rank;


  if (strcmp(flag, "-R") == 0) {
    for (i = 0; i < l_n; i++) {
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

int retreiveDimension(char* fileName){
  FILE* fp;
  char tmp[64];
  int rtn = 0;
  
  fp = fopen(fileName, "r");
  if (fp == NULL){
    printf("Emergency Abort: File Not Found!\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
  } 
  fgets(tmp, sizeof(tmp), fp);
  rtn = atoi(tmp);
  fclose(fp);

  return rtn;
}

void matrixMultiply(int local_N, int total_N, int M, double** mspace, double** ref, int maxR, double* global, int rank, int ncpu){
  int i, j, k, passes;
  int myNext, myLast, origin;
  int target[2] = {0, 0};

  myNext = rank + 1;
  myLast = rank - 1;
  origin = rank;

  if (rank == ncpu - 1) myNext = 0;
  if (rank == 0) myLast = ncpu - 1;
  int size = (maxR * M) + 1;
  printf("MaxR: %d M: %d\n",maxR, M);
  double* tempRef =  (double*)malloc(size * sizeof(double*));
  int tempSize = 0;
  
  tempRef[0] = (double)origin;
  
  for (i = 0; i < local_N; i++){
    for (j = 0; j < M; j++){
      tempRef[tempSize+1] = ref[i][j];
      tempSize++;
    }
  }
  
  for ( passes = 0; passes < ncpu; passes++){

    for (k = 0; k < M; k++){
      target[1] = tempRef[0];
      target[0] = 0;

      for (i = 0; i < local_N; i++){
    	for (j = 0; j < total_N; j++){
	  if (i == target[0] && j == target[1]){
	    mspace[i][total_N + k] += mspace[i][j] * tempRef[(M * i + k) + 1];
	    target[0]++;
	    target[1] = tempRef[0] + (i+1) * ncpu;
	  }
    	}
      }
    }

#ifdef DEBUG_PRINT_RECV
    MPI_Barrier(MPI_COMM_WORLD);

    for (i = 0; i < ncpu; i++){
    
      MPI_Barrier(MPI_COMM_WORLD);
    
      if( rank == i) {
	printf("Rank: %d :: ", rank);
	printVec(tempRef, size);
      }
      
      MPI_Barrier(MPI_COMM_WORLD);
    }
   
  
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // mspace[i][N+k] += mspace[i][j] * ref[M * i + k];
    MPI_Send(tempRef, size, MPI_DOUBLE, myNext, 0, MPI_COMM_WORLD);
    MPI_Recv(tempRef, size, MPI_DOUBLE, myLast, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    MPI_Barrier(MPI_COMM_WORLD);
    // printf("Moving on from pass: %d\n", passes);

  }
  
//   printf("Rank %d has the value of %d!!\n",inrank, myVal);
  
  free(tempRef);
}

/* void generateY(int N, int M, double** mspace, double** ref){ */
/*   int i, j, k; */
/*   double sum; */
  
/*   for (k = 0; k < M; k++){ */
/*     for (i = 0; i < N; i++){ */
/*       sum = 0; */
/*       for (j = 0; j < N; j++){ */
/*   	sum += mspace[i][j] * ref[i][k]; */
/*       } */
/*       mspace[i][N+k] = sum; */
/*     } */
/*   } */
/* } */

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
  int NxM[2] = {0,0};

  matFlag = argv[1];
  fileNameMat = argv[2];

  if (strcmp(matFlag, "-f") == 0){
    NxM[0] = retreiveDimension(fileNameMat);
  }else{
    NxM[0] = atoi(argv[2]);
  }
  
  refFlag = argv[3];
  fileNameRef = argv[4];
    
  if (strcmp(refFlag, "-F") == 0){
    NxM[1] = retreiveDimension(fileNameRef);
  }else{
    NxM[1] = atoi(argv[4]);
  }

#ifdef DEBUG_INPUTS
  printf("Inputs:  Flag1: %s n: %d fileMat: %s Flag2: %s m: %d fileRef: %s\n", matFlag, NxM[0], fileNameMat, refFlag,  NxM[1], fileNameRef);
#endif
  
  // Calculates which processor has what part of the global elements
  int local_n, tmp_n, biggest_n, R, i, j;

  double all_n[ncpu];
  
  for (i = 0; i < ncpu; i++){
    tmp_n= (int)( NxM[0]/ncpu);   // All processors have at least this many rows
    R = ( NxM[0]%ncpu);       // Remainder of above calculation
    if (i<R) {
      ++tmp_n;             // First R processors have one more element
    }
    all_n[i] = tmp_n;
  }

  local_n = all_n[rank];
  
  if (NxM[0] == ncpu){
    biggest_n = ((int)(NxM[0]/ncpu));
  }else{
    biggest_n = ((int)(NxM[0]/ncpu)) + 1;
  }
  
  double** A = (double**)malloc(local_n * sizeof(double*));
  double** r = (double**)malloc(local_n * sizeof(double*));
  double** x = (double**)malloc(local_n * sizeof(double*));
  
  for (i = 0; i < local_n; i++){
    A[i] = (double *)malloc(sizeof(double) * (NxM[0] + NxM[1]));
    r[i] = (double *)malloc(sizeof(double) * NxM[1]);
    x[i] = (double *)malloc(sizeof(double) * NxM[1]);
    assert(A[i] != NULL);
    assert(r[i] != NULL);
    assert(x[i] != NULL);
  }


  for (i = 0; i < local_n; i++) {
    for (j = 0; j < NxM[1]; j++){
      A[i][NxM[0] + j] = 0;
    }
  }
  

  createMatrix(local_n, NxM[0], A, matFlag, fileNameMat, rank, ncpu);
  createRefMatrix(local_n, NxM[0], NxM[1], r, refFlag, fileNameRef, rank, ncpu);
  matrixMultiply(local_n, NxM[0], NxM[1], A, r, biggest_n, all_n, rank, ncpu);
 
  /* for (k = 0; k < M; k++){ */
  /*   for (i = 0; i < N; i++){ */
  /*     sum = 0; */
  /*     for (j = 0; j < N; j++){ */
  /* 	sum += mspace[i][j] * ref[i][k]; */
  /*     } */
  /*     mspace[i][N+k] = sum; */
  /*   } */
  /* } */

  #ifdef DEBUG_PRINT_MAT
  
  MPI_Barrier(MPI_COMM_WORLD);

  for (i = 0; i < ncpu; i++){
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if( rank == i) {
      printf("Rank %d's matrix:\n", rank);
      printMat(A, local_n, NxM[0], NxM[1]);
      printf("Reference Matrix: \n");
      printMat(r, local_n, NxM[1], 0);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
   
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  #endif
  
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
