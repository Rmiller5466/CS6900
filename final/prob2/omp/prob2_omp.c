// cs6900
// Final - prob2_omp
// Ryan Miller
// w051rem
// Due 27 April 2021
// System = bender
// Compiler syntax = ./prob2_omp.compile prob2_omp
// Job Control File = prob2_omp.batch
// Additional File  = N/A
// Results file     = prob2_omp.txt

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>

#define slow 1000
#define eps 0.00000001

int main ( int argc, char **argv );
double local_minimum(double xguess, double yguess, double *x, double*y);
double Fxy(double x, double y);
double Dx(double x, double y, double dx);
double Dy(double x, double y, double dy);

// Global values
double r1=-8.1, r2=-2.1, r3=3, r4=7.2;
double r5=-7.9, r6=-2.5, r7=3.4, r8=6.8;
double X3, X2, X1, X0;
double Y3, Y2, Y1, Y0;
double xmin, xmax, ymin, ymax;

int main ( int argc, char **argv )
{
  // Global values
  X3=(r1 + r2 + r3 + r4);
  X2=(r1*r2 + r1*r3 + r1*r4 + r2*r3 + r2*r4 + r3*r4);
  X1=(r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4);
  X0=r1*r2*r3*r4;
  Y3=(r5 + r6 + r7 + r8);
  Y2=(r5*r6 + r5*r7 + r5*r8 + r6*r7 + r6*r8 + r7*r8);
  Y1=(r5*r6*r7 + r5*r6*r8 + r5*r7*r8 + r6*r7*r8);
  Y0=r5*r6*r7*r8;

  // Range
  xmax=10, xmin=-xmax;
  ymax=10, ymin=-ymax;

  // Initial dense grid to search over
  int nx=1000, ny=1000;
  double xstep=(xmax-xmin)/(nx-1);
  double ystep=(ymax-ymin)/(ny-1);

  int ix, iy;
  double xguess, yguess, xbest, ybest, x, y;
  double flocal, fmin=1e9; // initialize at a big value

  struct timeval start, end;
  unsigned long secs_used,micros_used;


  // Parallelize the below loop
  // brute force
  gettimeofday(&start, NULL);
#pragma omp parallel for shared(fmin, xbest, ybest) private(ix, iy, xguess, yguess, flocal)
  for (ix=0; ix<nx; ++ix){
    xguess=xmin+xstep*ix;
    for (iy=0; iy<ny; ++iy){
      yguess=ymin+ystep*iy;
      flocal=Fxy(xguess,yguess);
      if (flocal<fmin){
        fmin=flocal;
        xbest=xguess;
        ybest=yguess;
      }
    }
  } 
  gettimeofday(&end, NULL);

  
  
  printf("Brute Force answer \n");
  printf("x=%f y=%f min=%f \n",xbest,ybest,fmin);

  secs_used=(end.tv_sec - start.tv_sec); //avoid overflow by subtracting first
  micros_used= ((secs_used*1000000) + end.tv_usec) - (start.tv_usec);
  printf("Time for Brute Force answer %f sec \n\n",micros_used/1000000.0);


  // Parallelize the below loop
  // Initial grid to guess
  nx=4, ny=4;
  xstep=(xmax-xmin)/(nx-1);
  ystep=(ymax-ymin)/(ny-1);
  fmin=1e9; // initialize at a big value
  gettimeofday(&start, NULL);
#pragma omp parallel for shared(fmin, xbest, ybest) private(ix, iy, xguess, yguess, x, y, flocal)
  for (ix=0; ix<nx; ++ix){
    xguess=xmin+xstep*ix;
    for (iy=0; iy<ny; ++iy){
      yguess=ymin+ystep*iy;
      flocal=local_minimum(xguess,yguess,&x,&y);
      if (flocal<fmin){
        fmin=flocal;
        xbest=x;
        ybest=y;
      }
    }
  }
  gettimeofday(&end, NULL);
  printf("Smallest value %f\n",fmin);
  printf("Found at x=%f y=%f\n",xbest,ybest);
  printf("Within the range of x:%.2f %.2f\n",xmin,xmax);
  printf("                    y:%.2f %.2f\n",ymin,ymax);
  secs_used=(end.tv_sec - start.tv_sec); //avoid overflow by subtracting first
  micros_used= ((secs_used*1000000) + end.tv_usec) - (start.tv_usec);
  printf("Time for Steepest Descent Method %f sec \n",micros_used/1000000.0);

  return 0;
  //exit(0);
}

double local_minimum(double xguess, double yguess, double *x, double *y){
  double xstep=0.1;
  double ystep=0.1;
  
  double dx=Dx(xguess,yguess,xstep);
  double dy=Dy(xguess,yguess,ystep);
  double xnew=xguess-dx*xstep;
  double ynew=yguess-dy*ystep;
   if(xnew<xmin)
    xnew=xmin;
  if(xnew>xmax)
    xnew=xmax;
  if(ynew<ymin)
    ynew=ymin;
  if(ynew>ymax)
  ynew=ymax;

  double dxnew=Dx(xnew,ynew,xstep);
  double dynew=Dy(xnew,ynew,ystep);
  int count=0, itime=1;

  while (fabs(xguess-xnew)>eps && fabs(yguess-ynew)>eps && count<3){

    // If slope changed we pased the local min
    if ((dx*dxnew)<0){ // slope changed sign
      xstep=xstep/2.0;
    }else{
      xguess=xnew; 
    }
    if ((dy*dynew)<0){ // slope changed sign
      ystep=ystep/2.0;
    }else{
      yguess=ynew;
    }
    dx=Dx(xguess,yguess,xstep);
    dy=Dy(xguess,yguess,ystep);
    xnew=xguess-dx*xstep;
    ynew=yguess-dy*ystep;
    if(xnew<xmin)
      xnew=xmin;
    if(xnew>xmax)
      xnew=xmax;
    if(ynew<ymin)
      ynew=ymin;
    if(ynew>ymax)
      ynew=ymax;
    dxnew=Dx(xnew,ynew,xstep);
    dynew=Dy(xnew,ynew,ystep);
    ++itime;
    if(fabs(xguess-xnew)>eps || fabs(yguess-ynew)>eps){
      count=0;
    }else{
      ++count;
      xstep=xstep/2.0;
      ystep=ystep/2.0;
    }
  }
  //  printf("times %d \n",itime);
  // return location
  (*x)=xguess;
  (*y)=yguess;

  // return min vale at this location
  return(Fxy(xguess,yguess));
}

double Fxy(double x, double y){
  double fx, fy;
  usleep(slow);
  if (fabs(x)>xmax)
    return(0.0);
  if (fabs(y)>ymax)
    return(0.0);
  fx=cos(M_PI_2*x/xmax)*(x-r1)*(x-r2)*(x-r3)*(x-r4);
  fy=cos(M_PI_2*y/ymax)*(y-r5)*(y-r6)*(y-r7)*(y-r8);
  return(-fx*fy);
}

double Dx(double x, double y, double dx){
  return(((Fxy(x+dx/2.0,y)>Fxy(x-dx/2.0,y)) ? 1.0 : -1.0 ));
}

double Dy(double x, double y, double dy){
  return(((Fxy(x,y+dy/2.0)>Fxy(x,y-dy/2.0)) ? 1.0 : -1.0 ));
}

