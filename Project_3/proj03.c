// cs6900
// Project 03
// Ryan Miller
// File Template Provided by: Dr. John Nehrbass
// w051rem
// Due 12 Feb 2021
// System = bender
// Compiler syntax = ./compile.sh proj03
// Job Control File = proj03.sbatch
// Additional File  = NA
// Results file     = proj03.txt

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#define DEBUG

// DO NOT CHANGE THIS FUNCTION!!!!!
// function to return a random variable between [0,1)
double randVar(){

  // Notes:
  // 2147483647 is the largest possible value of rand()
  // printf("max=%llu\n",RAND_MAX);

  return (double)((rand()%2000000000)/2000000000.0);
}


/* function to return M and P*/
unsigned long* calMandP(int N ) {

  //Save the results in an integer array
  static unsigned long MnP[2];

  // Constant for small cube
  // 1/sqrt(3)
  double rsq3=0.57735026919;

  // initialize counters each time the function is called
  unsigned long M=0, P=0;
  
  // random variables
  double x, y, z, rs;
  int randomPnts;

   // Pick N random (x,y) points
  for(randomPnts = 0; randomPnts < N; ++randomPnts) {

    x=randVar();
    y=randVar();
    z=randVar();

    // radius squared
    rs=x*x+y*y+z*z;

#ifdef DEBUG
    // test random values in function
    printf("%d x=%f y=%f z=%f rs=%f\n",randomPnts,x,y,z,rs);
#endif 

    if (rs<1){                  
      // Inside circle
      ++M;
      if (x<rsq3 && y<rsq3 && z<rsq3){ 
	// inside inner square
	++P;
      }
    }
  }

  MnP[0]=M;  
  MnP[1]=P;  
  return MnP;
}


int main (int argc, char *argv[]){
  int start_time = time(NULL);

  MPI_Init(NULL, NULL);

  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  char processor_name[MPI_MAX_PROCESSOR_NAME];

  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  
  srand(time(NULL)*(double)world_rank);
  // For convenience I passed "N" on the command line
  // Not checking for valid input
  char *a = argv[1];
  char *b = argv[2];
  int N = atoi(a);
  double E = atof(b);
  int running = 1;
  // The Max value of unsigned long is 18,446,744,073,709,551,615
  unsigned long* MnP;
  unsigned long MnP_temp[2] = {0};
  
  if (world_rank == 0){
    
    unsigned long MnP_total[2] = {0};
    double sqrt3 = 1.73205080757;
    double epsilon = E;
    int n_modifier = world_size - 1;
    int MAX_RUNS = 100000;
    int current_runs=1;
    double cube2sphere = 0, sphere2cube = 0, est_value1 = 0, est_value2 = 0;
    int finished_1 = 0, finished_2 = 0, identifier = 0;
    MnP = calMandP(0);

    while ((current_runs < MAX_RUNS) && running == 1){
      
      MPI_Reduce(MnP,MnP_temp,2,MPI_UNSIGNED_LONG,MPI_SUM,0,MPI_COMM_WORLD);
      MnP_total[0] += MnP_temp[0];
      MnP_total[1] += MnP_temp[1];
      
      if (finished_1 == 0){
	
	cube2sphere=6.0*((double)MnP_total[0]/(double)(n_modifier * current_runs * N));
	if (fabs(cube2sphere - est_value1) > epsilon) {
	    est_value1 = cube2sphere;
	  }else{
	  finished_1  = 1;
	  identifier = 1;
	}
      }
      

      if (finished_2 == 0){
	
	sphere2cube=(2.0*(((double)MnP_total[0]/(double)MnP_total[1])/sqrt3));
	if ( fabs(sphere2cube - est_value2) > epsilon) {
	    est_value2 = sphere2cube;
	  }else{
	  finished_2  = 1;
	  identifier = 2;
	}
      }
      
      if (finished_1 == 1 || finished_2 == 1) {
	running = 0;
      }
     
      current_runs += 1;
      if (current_runs >= MAX_RUNS) {
	running = 0;
      }
      
      MPI_Bcast(&running, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    int end_time = time(NULL);
    int total_time = end_time - start_time;
    
    
    printf("\nProgram Information:\nTotal Processors: %d Total Time: %d N: %d\n\n",world_size, total_time, N);
    printf("Epsilon: %s\nIterations: %d\nM = %d P = %d\nEstimation of Pi (Cube:Sphere): %f\nEstimation of Pi (Sphere:Cube): %f\nMethod %d converged\n",
   	   b, current_runs, MnP_total[0], MnP_total[1], est_value1, est_value2, identifier);
    
    
  }else{
   
    while( running == 1){
        
      MnP=calMandP(N);
      
      
      MPI_Reduce(MnP,MnP_temp,2,MPI_UNSIGNED_LONG,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&running, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
  }
  MPI_Finalize();
}
