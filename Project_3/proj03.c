//
// cs4900
// Project 03
// Dr. John Nehrbass
// w006jwn
// Due 12 Feb 2021
// System = bender
// Compiler syntax = ./compile.sh proj03
// Job Control File = proj03.sbatch
// Additional File  = NA
// Results file     = proj03.txt
//

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

    if (rs<1){                   // or should it be rs<=1
      // Inside circle
      ++M;
      if (x<rsq3 && y<rsq3 && z<rsq3){     // or should it be <=
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
  MPI_Init(NULL, NULL);

  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  char processor_name[MPI_MAX_PROCESSOR_NAME];

  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);



  // Use current time as seed for random generator 
  srand(time(NULL)*(double)world_rank); 
  
  // For convenience I passed "N" on the command line
  // Not checking for valid input
  char *a = argv[1];
  int N = atoi(a);
  
  // The Max value of unsigned long is 18,446,744,073,709,551,615
  unsigned long* MnP;   // Store in vector M=MnP[0] and P=MnP[1]
  
  // Return M and P for this case
  MnP=calMandP(N);
  
  // Results
   
  double sqrt3=1.73205080757;
  double p1=6.0*((double)MnP[0]/(double)N);
  // double p2 = 3.0*((double)MnP[0]/(double)MnP[1]);
  double p2=(2.0*(((double)MnP[0]/(double)MnP[1])/sqrt3));
  // double p1=(2.0*((double)MnP[0]/(double)N))/sqrt3;
	     printf("Reporting From: Processor %s, Rank %d Of %d\nM=%llu P=%llu\nPi1: %f Pi2: %f\n",  processor_name, world_rank, world_size, MnP[0], MnP[1],p1,p2);
  // printf(,p1,p2);

  MPI_Finalize();// Finalize the MPI environment.
}

// Get the number of processes
// Get the rank of the process

// Get the name of the processor 

//printf("Hello world from processor %s, rank %d out of %d processors\n",  processor_name, world_rank,  world_size);



