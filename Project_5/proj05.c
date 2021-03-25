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

/*
  Debug segments used for troubleshooting, viewing inputs and matrices
  Uncomment each #define desired
*/

#define DEBUG_PRINT_MAT
//#define DEBUG_PRINT_REF
//#define DEBUG_PRINT_RECV
//#define DEBUG_INPUTS

/*
  Method: printVec
  Inputs: (double*) vector - pointer to vector
          (int) N - Length of vector
 
  Functionality: Print to screen the passed Vector
*/

void printVec(double* vector, int N){
  int i;

  for (i = 0; i < N; i++) {
    if (i == (N - 1)){
      printf("%.20f\n", vector[i]);
    }else{
      printf("%.20f ", vector[i]);
    }
  }
}

/*
  Method: printMat
  Inputs: (double**) mat - pointer to a 2 dimensional matrix
          (int) local_N - Number of rows
	  (int) total_N - Number of columns (Also the global N)
	  (int) M - Number of solution columns
  
  Functionality: Print to screen the passed Matrix
*/

void printMat(double** mat, int local_N, int total_N, int M){
  int i,j;
  
  for (i = 0; i < local_N; i++) {
    for (j = 0; j < total_N+M; j++) {
      if (j == total_N){
	// Divider between solution columns and x columns
	printf(": ");
      }
      if (j == (total_N+M - 1)){
  	printf("%.3f\n", mat[i][j]);
      }else{
  	printf("%.3f ", mat[i][j]);
      }
    }
  }
}


/*
  Method: createMatrix
  Inputs: (int) local_N - Number of rows
	  (int) total_N - Number of columns (Also the global N) 
          (double**) mspace - Stands for matrix space, where the data generated will be stored
          (char*) flag - Character string representing how to create the matrix
	  (char*) filename - If reading from a file, this variable will hold the name of the file
	  (int) rank - Processor number
	  (int) ncpu - Total number of processors
  
  Functionality: Populates an allocated matrix space with values depending on chosen 'mode' of 
                 creation. The flag will be -r, -d, or -f.
		 -r : Initializes the matrix with random values
		 -d : Initializes the matrix as a diagonal matrix
		 -f : Initializes the matrix from a provided file
*/

void createMatrix(int local_N, int total_N, double** mspace, char* flag, char* filename, int rank, int ncpu){ 

  int i, j; 
  int offset = rank;

  //First if branch | Creates a matrix of random nubmers
  if (strcmp(flag, "-r") == 0) {
    for (i = 0; i < local_N; i++) {
      for (j = 0; j < total_N; j++) {
	mspace[i][j] = (double)((rand() % 100) + rank);  
      }
    }

  // Second if branch | Creates diagonal based matrix using an offset depending on the rank
  }else if (strcmp(flag, "-d") == 0){
    for (i = 0; i < local_N; i++) {

      offset = rank + (i * ncpu);

      for (j = 0; j < total_N; j++) {

	if ( offset == j) {
	  mspace[i][j] =  offset + 1;
	}else{
	  mspace[i][j] = 0;
	}
      }
    }

  // Final if branch | Reads in matrix from a given file
  }else if(strcmp(flag, "-f") == 0) {
    
    FILE* fp;
    
    fp = fopen(filename, "r"); 

    // Filename was invalid... exit program
    if (fp == NULL){
      printf("[createMatrix] Emergency Abort: File Not Found!\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(1);
    }

    // Hopefully the input string size is no bigger than 4096 
    char st[4096];
    char* num;
    char* eof;

    // Skips the first line (We already got N in main() )
    fgets(st, sizeof(st), fp);
    eof = fgets(st,sizeof(st), fp);
 
    // If there was not a second line, exit program
    if (eof == NULL) {
      printf("[createMatrix] Emergency Abort: Incorrect File Format!\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(1);
    }
    
    fclose(fp);
    // Print the received Matrix String
    //printf("String: %s\n",st);
    
    int counter = 0;
      
    num = strtok(st, " ");
    for (i = 0; i < total_N; i++) { 
      for (j = 0; j < total_N; j++) {
	
	// num should never be NULL, if it is, exit program
	if (num == NULL){
	  printf("[createMatrix] Emergency Abort: Incorrect File Format!\n");
	  MPI_Abort(MPI_COMM_WORLD, 1);
	  exit(1);
	}
	if (i == offset){
	  mspace[counter][j] = atof(num);	
	}
	
	//Updates the offset so that the rank can read the correct row
	if (i == offset && j == total_N - 1){
	  counter++;
	  offset = rank + (counter * ncpu);
	}
	num = strtok(NULL, " ");
      }
    }

  // If an invalid flag was passed in, exit program
  }else{
    printf("[createMatrix] Emergency Abort: Invalid Input Flag (Must use -r, -d, or -f!)\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
  }
}


/*
  Method: createRefMatrix
  Inputs: (int) local_N - Number of rows
	  (int) total_N - Number of columns (Also the global N)
	  (int) M - Number of reference / solution columns
          (double**) ref - Memory space where the data generated will be stored, the reference matrix
          (char*) flag - Character string representing how to create the matrix
	  (char*) filename - If reading from a file, this variable will hold the name of the file
	  (int) rank - Processor number
	  (int) ncpu - Total number of processors
  
  Functionality: Populates an allocated matrix space with values depending on chosen 'mode' of 
                 creation. The flag will be -R or -F.
		 -R : Initializes the reference matrix with random values
		 -F : Initializes the reference matrix from a provided file
*/

void createRefMatrix(int local_N, int total_N, int M,  double** ref, char* flag, char* filename, int rank, int ncpu){ 
  int i, j;
  int offset = rank;

  // Creates a reference matrix of random numbers
  if (strcmp(flag, "-R") == 0) {
    for (i = 0; i < local_N; i++) {
      for (j = 0; j < M; j++) {
	ref[i][j] = (double)((rand() % 50) + rank);  
      }
    }

  // Creates the reference matrix from the provided input file
  }else if(strcmp(flag, "-F") == 0) {
    
    FILE* fp;
    
    fp = fopen(filename, "r"); 

    // If the file was not found, exit the program
    if (fp == NULL){
      printf("[createRefMatrix] Emergency Abort: File Not Found!\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(1);
    }

    // Again, we hope the input string sizes are not greater than 4096
    char st[4096];
    char* eof;
    char* num;
    int counter = 0;

    fgets(st, sizeof(st), fp);
    for (i = 0; i < M; i++){
      eof = fgets(st,sizeof(st), fp);
      
      // If the file format is incorrect, exit the program
      if (eof == NULL) {
	printf("[createRefMatrix] Emergency Abort: Input File Line Count Didn't Match M!\n");
	MPI_Abort(MPI_COMM_WORLD, 1);
	exit(1);
      }
      num = strtok(st, " ");
      
      for (j = 0; j < total_N; j++) {
	if (j == offset){
	  ref[counter][i] = atof(num);
	  counter++;
	  offset = rank + (counter * ncpu);
	}
	num = strtok(NULL, " ");
      }
      counter = 0;
      offset = rank;
    }
    fclose(fp);

  // If the input flag given is invalid, exit the program
  }else{
    printf("[createRefMatrix] Emergency Abort: Invalid Input Flag (Must use -R, or -F!)\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
  }
  
}



/*
  Method: retreiveDimension
  Inputs: (char*) filename - Name of the file to be read from

  Functionality: Reads from the provided file and returns a single value on the first line of the document  
*/

int retreiveDimension(char* filename){
  FILE* fp;
  char tmp[64];
  int rtn = 0;

  fp = fopen(filename, "r");
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

/*
  Method: matrixMultiply
  Inputs: (int) local_N - Number of rows
	  (int) total_N - Number of columns (Also the global N) 
          (int) M - Number of reference / solution columns
	  (double**) mspace - The 2d vector representing a rank's A matrix
	  (double**) ref - The 2d vector representing a ranks' reference matrix
	  (int) maxR - Represents the maximum amount of rows any rank will have
	  (int) rank - Processor number
	  (int) ncpu - Total number of processors
  
  Functionality: Populates an allocated matrix space with values depending on chosen 'mode' of 
                 creation. The flag will be -r, -d, or -f.
		 -r : Initializes the matrix with random values
		 -d : Initializes the matrix as a diagonal matrix
		 -f : Initializes the matrix from a provided file
*/

void matrixMultiply(int local_N, int total_N, int M, double** mspace, double** ref, int maxR, int rank, int ncpu){
  int i, j, k, passes;
  int counter;
  int myNext, myLast, origin;
  
  //  printf("Local: %d, Total: %d, M: %d maxR: %d\n", local_N, total_N, M, maxR);


  myNext = rank + 1;
  myLast = rank - 1;
  origin = rank;

  if (rank == ncpu - 1) myNext = 0;
  if (rank == 0) myLast = ncpu - 1;
  
  int size = (maxR * M) + 1;
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

      for (i = 0; i < local_N; i++){
	counter = 0;
      	for (j = 0; j < total_N; j++){
      	  if (j == tempRef[0] + counter * ncpu){
      	    mspace[i][total_N + k] += mspace[i][j] * tempRef[(M * counter + k) + 1];
	    counter++;
	  }
      	}
      } 
    }

/* #ifdef DEBUG_PRINT_RECV */
/*     MPI_Barrier(MPI_COMM_WORLD); */

/*     for (i = 0; i < ncpu; i++){ */
    
/*       MPI_Barrier(MPI_COMM_WORLD); */
    
/*       if( rank == i) { */
/* 	printf("Rank: %d :: ", rank); */
/* 	printVec(tempRef, size); */
/*       } */
      
/*       MPI_Barrier(MPI_COMM_WORLD); */
/*     } */
   
  
/*     MPI_Barrier(MPI_COMM_WORLD); */
/* #endif */

    MPI_Send(tempRef, size, MPI_DOUBLE, myNext, 0, MPI_COMM_WORLD);
    MPI_Recv(tempRef, size, MPI_DOUBLE, myLast, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
  }
  
  free(tempRef);
}

void pivot(double** mat, int local_N, int total_N, int M, int currentCol, int rank, int ncpu){
  int i;
  int max_Local_Index, maxRank;
  int swapRank, swap_Local_Index;
  int size = total_N + M;

  struct { 
    double value; 
    int   index; 
  } localMax, globalMax; 


  // Local Max
  localMax.value = 0; 
  localMax.index = 0; 
  for (i=0; i < local_N; i++){ 
    if (rank + ncpu * i >= currentCol){
      if (fabs(localMax.value) < fabs(mat[i][currentCol])) { 
	localMax.value = mat[i][currentCol]; 
	localMax.index = i; 
      } 
    }
  }
  
  localMax.index = rank + ncpu * localMax.index; 

  MPI_Reduce(&localMax, &globalMax, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD );
  MPI_Bcast(&globalMax, 1, MPI_DOUBLE_INT, 0, MPI_COMM_WORLD);

  //All Ranks now have the global max

  maxRank = globalMax.index % ncpu;
  max_Local_Index = globalMax.index / ncpu;

  swapRank = currentCol % ncpu;
  swap_Local_Index = currentCol / ncpu;

  
  if (swapRank != maxRank){
    
    double* tempRow;
    
    if (rank == swapRank){
      
      tempRow = mat[swap_Local_Index];
      
      MPI_Send(tempRow, size, MPI_DOUBLE, maxRank, 0, MPI_COMM_WORLD);
      MPI_Recv(tempRow, size, MPI_DOUBLE, maxRank, 0, MPI_COMM_WORLD,
    	       MPI_STATUS_IGNORE);
      
      mat[swap_Local_Index] = tempRow;
    }

    if (rank == maxRank){

      tempRow = mat[max_Local_Index];
      MPI_Send(tempRow, size, MPI_DOUBLE, swapRank, 0, MPI_COMM_WORLD);
      MPI_Recv(tempRow, size, MPI_DOUBLE, swapRank, 0, MPI_COMM_WORLD,
    	       MPI_STATUS_IGNORE);
      mat[max_Local_Index] = tempRow;
    }
  }else if (swapRank == maxRank && rank == maxRank && max_Local_Index != swap_Local_Index){

     double* tempRow;
     tempRow = mat[swap_Local_Index];
     mat[swap_Local_Index] = mat[max_Local_Index];
     mat[max_Local_Index] = tempRow;
   

  } 
  /* if (rank == 0) {  */
  /*   printf("Pivot Col: %d, Maxval=%f maxrank=%d maxindexG=%d maxindexL=%d\n",currentCol,maxVal,maxRank,max_Global_Index, max_Local_Index); */
  /*   printf("Pivot Col: %d, currank=%d curindexG=%d curindexL=%d\n",currentCol,currentRank,current_Global_Index, current_Local_Index); */

  /* } */
}

void GE(int local_N, int total_N, int M, double** mat,  int rank, int ncpu){
  int i, base, row, col;
  double cancelVal;
  int baseRank, base_Local_Index;
  int size = total_N + M;
  double tempR[size];

  
  for (base = 0; base < total_N-1; base++) {
    pivot(mat, local_N, total_N, M, base, rank, ncpu);

    baseRank=base%ncpu;
    base_Local_Index = (base-baseRank)/ncpu;

    if (rank == baseRank){
      for (i = 0; i < size; i++){
      tempR[i] = mat[base_Local_Index][i];
      }
    }
    MPI_Bcast(&tempR, size, MPI_DOUBLE, baseRank, MPI_COMM_WORLD);	
    

    for (row = 0; row < local_N; row++){
      if (row * ncpu + rank > base){
	cancelVal = mat[row][base]/tempR[base];
	for (col = 0; col < size; col++){
	  mat[row][col] = mat[row][col] + (-1 * cancelVal * tempR[col]);
	}
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

}

void BS(int local_N, int total_N, int M, double** mat, double** xVec,  int rank, int ncpu){
  int base, row, i;
  int baseRank, base_Local_Index;
  double tempX[M];

  for ( base = total_N-1; base >= 0; base--){
    baseRank=base%ncpu;
    base_Local_Index = (base-baseRank)/ncpu;
    
    if (rank == baseRank){
      for (i = 0; i < M; i++){
	//	printf("[%d] tempX: %d will be: %f divided by: %f\n",rank, i, mat[base_Local_Index][total_N+i],  mat[base_Local_Index][base]);
	tempX[i] = mat[base_Local_Index][total_N+i] / mat[base_Local_Index][base];
	
	//printf("[%d] tempX: %d is now: %f on base: %d\n",rank, i, tempX[i], base);

	xVec[base_Local_Index][i] = tempX[i];
	      
      }
    }
    
    MPI_Bcast(&tempX, M, MPI_DOUBLE, baseRank, MPI_COMM_WORLD);	
    
    for (row = 0; row < local_N; row++){
      if (row * ncpu + rank < base){
    	for (i = 0; i < M; i++){
	  // printf("[%d] doing %f - %f * %f\n",rank, mat[row][total_N+i], tempX[i], mat[row][base]);
  	  mat[row][total_N+i] = mat[row][total_N+i] - tempX[i] * mat[row][base];
    	}
	 mat[row][base] = 0;
      }
    }
    
  }
}


int main(int argc, char** argv) {

  int rank, ncpu;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ncpu);

  srand(time(NULL)*(double)(rank+1));

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
  double* errorVec = (double*)malloc(NxM[1] * sizeof(double*));
  double* epsilon = (double*)malloc(NxM[1] * sizeof(double*));

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

  for(i = 0; i < NxM[1]; i++){
    errorVec[i] = 0;
    epsilon[i] = 0;
  }
  

  createMatrix(local_n, NxM[0], A, matFlag, fileNameMat, rank, ncpu);
  createRefMatrix(local_n, NxM[0], NxM[1], r, refFlag, fileNameRef, rank, ncpu);
  matrixMultiply(local_n, NxM[0], NxM[1], A, r, biggest_n, rank, ncpu);

  GE(local_n, NxM[0], NxM[1],  A, rank, ncpu);
  BS(local_n, NxM[0], NxM[1], A, x, rank, ncpu);

  for( i = 0; i < local_n; i++){
    for (j = 0; j < NxM[1]; j++){
      errorVec[j] += pow(x[i][j] - r[i][j], 2);
    }
  }

  MPI_Reduce(errorVec,epsilon,NxM[1],MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  if (rank == 0){
    for( i = 0; i < NxM[1]; i++){
      epsilon[i] = sqrt(epsilon[i]) / NxM[0];
    }
    if (epsilon[0] > 1){
      printf("Error! System has no Solution!\n");
    }else{
      
      printVec(epsilon, NxM[1]);
    }
  }

/* #ifdef DEBUG_PRINT_MAT */
  
/*   MPI_Barrier(MPI_COMM_WORLD); */
/*   for (i = 0; i < ncpu; i++){ */
    
/*     MPI_Barrier(MPI_COMM_WORLD); */
    
/*     if( rank == i) { */
/*       printf("Rank %d's matrix:\n", rank); */
/*       printMat(A, local_n, NxM[0], NxM[1]); */

/*       printf("\nReference Matrix: \n"); */
/*       printMat(r, local_n, NxM[1], 0); */

/*       printf("Solution Matrix: \n"); */
/*       printMat(x, local_n, NxM[1], 0); */
      
/*       printf("ErrorVec:\n"); */
/*       printVec(errorVec, NxM[1]); */

/*       printf("Epsilon:\n"); */
/*       printVec(epsilon, NxM[1]); */
/*     } */

/*     MPI_Barrier(MPI_COMM_WORLD); */
/*   } */
   
  
/*   MPI_Barrier(MPI_COMM_WORLD); */
  
/* #endif */



  for (i = 0; i < local_n; i++) {
    free(A[i]);
    free(r[i]);
    free(x[i]);
  }
  free(A);
  free(r);
  free(x);
  free(errorVec);
  free(epsilon);
  

  MPI_Finalize(); 
}



