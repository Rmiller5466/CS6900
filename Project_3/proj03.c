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
  double start_time = time(NULL);

  MPI_Init(NULL, NULL);

  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  char processor_name[MPI_MAX_PROCESSOR_NAME];

  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
   
  // For convenience I passed "N" on the command line
  // Not checking for valid input
  char *a = argv[1];
  int N = atoi(a);
 
  // The Max value of unsigned long is 18,446,744,073,709,551,615
  unsigned long* MnP,MnP_temp,MnP_total;   // Store in vector M=MnP[0] and P=MnP[1]

  if (world_rank == 0){
    
    double sqrt3=1.73205080757;
    double epsilon = 0.001;
    int MAX_RUNS = 1000;
    int current_runs = 0;
    double cube2sphere = 0, sphere2cube = 0, est_value1 = 0, est_value2 = 0;
    int finished_1 = 0, finished_2 = 0;

    while (current_runs < MAX_RUNS){
      
      if (finished_1 == 1 && finished_2 == 1) {
	//TODO: Signal to stop
      }
      
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Reduce(MnP,MnP_temp,2,MPI_UNSIGNED_LONG,MPI_SUM,0,MPI_COMM_WORLD);
      
      MnP_total[0] += MnP_temp[0];
      MnP_total[1] += MnP_temp[1];
      
      if (finished_1 == 0){
	
	cube2sphere=6.0*((double)MnP_total[0]/(double)N);
	if ((cube2sphere - est_value1) > epsilon) || (( est_value1 - cube2sphere) > epsilon) {
	    est_value1 = cube2sphere;
	  }else{
	  finished_1  = 1;
	}
      }
      

      if (finished_2 == 0){
	
	sphere2cube=(2.0*(((double)MnP_total[0]/(double)MnP_total[1])/sqrt3));
	if ((sphere2cube - est_value2) > epsilon) || (( est_value2 - sphere2cube) > epsilon) {
	    est_value2 = sphere2cube;
	  }else{
	  finished_2  = 1;
	}
      }
    }
    
    //TODO: Finalize and print Final Results


  }else{
    
    // Use current time as seed for random generator
    srand(time(NULL)*(double)world_rank);
    
    //TODO: Set up M and P and return to rank 0
    //  MnP=calMandP(N);
  }
  
  //USE LATER:
  //  printf("Reporting From: Processor %s, Rank %d Of %d\nM=%llu P=%llu\nPi1: %f Pi2: %f\n",  processor_name, world_rank, world_size, MnP[0], MnP[1],p1,p2);
  
}
