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

  // Constant for small square 
  // Reciprocal of sqrt(2) or 1/sqrt(2)
  double rsq2=0.70710678118654752440;

  // initialize counters each time the function is called
  unsigned long M=0, P=0;
  
  // random variables
  double x, y, rs;
  int randomPnts;

   // Pick N random (x,y) points
  for(randomPnts = 0; randomPnts < N; ++randomPnts) {

    x=randVar();
    y=randVar();

    // radius squared
    rs=x*x+y*y;

#ifdef DEBUG
    // test random values in function
    printf("%d x=%f y=%f rs=%f\n",randomPnts,x,y,rs);
#endif 

    if (rs<1){                   // or should it be rs<=1
      // Inside circle
      ++M;
      if (x<rsq2 && y<rsq2){     // or should it be <=
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

  // Use current time as seed for random generator 
  srand(time(0)); 
  
  // For convenience I passed "N" on the command line
  // Not checking for valid input
  char *a = argv[1];
  int N = atoi(a);
  
  // The Max value of unsigned long is 18,446,744,073,709,551,615
  unsigned long* MnP;   // Store in vector M=MnP[0] and P=MnP[1]
  
  // Return M and P for this case
  MnP=calMandP(N);
  
  // Results
  printf("M=%llu P=%llu\n",MnP[0],MnP[1]);
  
  return 0;
}




