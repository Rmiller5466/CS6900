# CS6900 - HPC & Parallel Programming
Location where projects being done throughout this class will be stored.
1. Project 1 - Warming up with Bash/Linux
   * Writing a basic bash script
   * Uses input passed with the program call
   * Calculates the sum and average of input integers
   * Includes -h and -v page
   * Validates inputs and handles errors correctly

2. Project 2 - Getting familiar with SLURM
   * Write a basic SLURM Job Control File
   * Utilizes Job Arrays
   * Decrypt program supplied by Dr. John Nehrbass
   * Int1 and Int2 are passed to Decrypt (along with student ID)
   * Int1 has 16 possiblites, Int2 has 512
   * All these combinations must be tested to see which codes access to data

3. Project 3 - Estimate Pi Using Parallel Programming
   * Working with a Cube inscribed in a Sphere which is also inscribed in another Cube
   * Radius of Sphere is 1
   * MPI program starts and estimates N points within the shapes (Positive numbers only)
   * These points are generated on P processors
   * M = Points in Sphere, P = Points in both Sphere and Inner Cube
   * Using M and P, an estimation of Pi is calculated on Rank 0
   * The program continues generating points and adding to the global total
   * Once Pi changes by less than epsilon, all Ranks terminate, results are reported
   * Research how message size impacts message passing between ranks

4. Project 4 - Parallel Dot Product And Matrix Multiplication
   * Inputs include vector lengths, offsets, and coefficients
   * Partition vectors across processors
   * Calculate both dot product and matrix resulting from matrix multiplication 
   * Rank calculates with only the data that it has been given
   * Run these with many different inputs and sizes, plot and document the results
   * Automate all 31 runs

