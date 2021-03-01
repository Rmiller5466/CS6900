#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <mpi.h>


double dot(double *v1, double *v2, int N) {
  int i;
  double sum=0.0;
  for (i = 0; i < N; i++) {
    double temp = v1[i]*v2[i];
    //printf("i value: %d, %f ",i,temp);
    sum += v1[i]*v2[i];
  }
  //printf("\n\n");
  return sum;
}

double* calcMatrix(double **matrix, int Nr, int l2, int offset){
  static double sums[2] = {0};
  int i,j;

  for (i = 0; i < Nr; i++){
    for (j = 0; j < l2; j++){
      //printf("%f ", matrix[i][j]);
      if (j == (offset-1+i)){
	//printf("%f\n", matrix[i][j]);
	sums[1] += matrix[i][j];
      }else{
	sums[0] += matrix[i][j];
      }
      
    }
    //printf("\n");
  }
  //  printf("Sums from offset: %d, 0(ND)=%f , 1(D)=%f", offset, sums[0], sums[1]);
  return sums;
}

 


int main(int argc, char** argv) {

  int rank, ncpu, name_len;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ncpu);
  


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
  
  int MAT_CONST[] = {0, 1000, 100000};
  int OFFSET_CONST[] = {22, 45, 75, 90, 135};
  int runs;
  
  int L1, L2, A, B, Ap, Bp;
  L2 = MAT_CONST[0];
  for (runs = 0; runs < 31; runs++){
    if (runs < 9){

      L1 = L1_CONST[runs % 9];
      A = 1; B = -1; Ap = OFFSET_CONST[1]; Bp = OFFSET_CONST[1];
    
    }
    if (runs > 8 && runs <18){
    
      L1 = L1_CONST[runs % 9];
      A = 1; B = (int)(L1/4); Ap = OFFSET_CONST[0]; Bp = OFFSET_CONST[3];

    }
    if (runs > 17 && runs < 27){
      
      L1 = L1_CONST[runs % 9];
      A = 1; B = (int)(2.2*L1); Ap = OFFSET_CONST[2]; Bp = OFFSET_CONST[4];
    
    }
    if (runs > 26 && runs < 29){
    
      L1 = MAT_CONST[runs % 2]; L2 = L1;
      A = 1; B = 1; Ap = OFFSET_CONST[4]; Bp = OFFSET_CONST[4];
    
    }
    if (runs > 28 && runs < 31){
    
      L1 = MAT_CONST[runs % 2]; L2 = L1;
      A = 1; B = (int)(2.4*L2); Ap = OFFSET_CONST[1]; Bp = OFFSET_CONST[4];
   
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
  
    //  printf("Pi: %f Dx: %f Dy: %f norm1: %f norm2: %f a: %f b: %f\n", pi, dx, dy, norm1, norm2, a, b);
  
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
  
    //(i = row-1; i < row+Nr-1; i++)
  
    for (i = 0; i < Nr; i++)  {
      Xi = dx/2.0 + (row-1+i)  * dx;
      v1[i] = norm1*cos((A*Xi)+a);
      v2[i] = norm1*sin((B*Xi)+b);
    }

    for (i = 0; i < L2; i++){
      Yi = dy/2.0 + i * dy;
      v3[i] = norm2*cos((B*Yi)+b);
    }
  
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
	// printf("v1 J-%d VALUE: %f\nv3 I-%d VALUE: %f\nMATRIX-I: %d J: %d VALUE:%f\n\n",j,v1[j],i,v3[i],i,j,Matrix[i][j]);
      }
    }
  

  
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

    /*    if( rank == 1) { */
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
   
    /*    if( rank == 2) { */
    /*     printf("Rank: %d, Matrix: \n", rank); */
    /*     for (i = 0; i < Nr; i++) { */
    /*       for ( j = 0; j < L2; j++){ */
    /* 	 printf("%f ", Matrix[i][j]); */
    /*       } */
    /*       printf("\n"); */
    /*     } */
    /*     printf("\n"); */
    /*   } */
  
    /* MPI_Barrier(MPI_COMM_WORLD); */
  

    double *localMatrixTotals = calcMatrix(Matrix, Nr, L2, row);

    double dotSumRecv;
    double matrixTotalsRecv[2] = {0};

    if (rank == 0){
      MPI_Reduce(&localDotSum,&dotSumRecv,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(localMatrixTotals,matrixTotalsRecv,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   
      printf("L1: %d L2:%d A: %d B: %d Ap: %d Bp: %d\n",L1, L2, A, B, Ap, Bp);
      printf("Dot Product of V1:V2 was: %f\n",dotSumRecv);
      printf("Matrix Total (Non-Diagonal): %f, Matrix Diagonal Total: %f\n\n",matrixTotalsRecv[0], matrixTotalsRecv[1]);
   
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
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}
