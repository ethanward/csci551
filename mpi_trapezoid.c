/*

  Author: Ethan Ward
  ID: 005704518

  This programs utilizes open_mpi to perform the Trapezoidal Rule and
  collect/analyze empirical timing data for the purpose of testing
  the effects of using multiple processes for a program.  

*/

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Approc Val: 4003.7209001513997464
// True Val: 4003.7209001513268265


/*
  Input: Any double value
  Ouput: The solution to the function cos(x/3) - 2*cos(x/5) + 5*sin(x/4) + 8
*/
double f(double x) { 
  return (cos(x/3) - 2*cos(x/5) + 5*sin(x/4) + 8); 
}

/*
  Input: The a, b, n, and h to the Trapizoidal Rule
  Output: The area under the function f(double x) curve using the trapezoid rule with
    the given a, b, n, and h values. 
*/
double Trap(double left, double right, int trap_count, double base) {
  double estimate, x;
  int i;

  estimate = (f(left) + f(right)) / 2.0;
  for(i = 1; i <= trap_count-1; i++) {
    x = left + i*base;
    estimate += f(x);
  }

  estimate = estimate*base;

  return estimate;
}

/*
  Input: The processor running the function, the number of processors, pointers 
    to ouput variables.
  Output: none. 

  Copied from Pacheco Program 3.5
*/
void Get_input(int my_rank, int comm_sz, double* a_p, double* b_p, int* n_p) {
   int dest;

   if (my_rank == 0) {
      printf("Enter a, b, and n\n");
      scanf("%lf %lf %d", a_p, b_p, n_p);
      for (dest = 1; dest < comm_sz; dest++) {
        MPI_Send(a_p, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
        MPI_Send(b_p, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
        MPI_Send(n_p, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
      } 
   } 
   else { /* my_rank != 0 */
      MPI_Recv(a_p, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);
      MPI_Recv(b_p, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);
      MPI_Recv(n_p, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);
   } 
} 

/*
  After getting user input for the trapezoidal rule, this function calculates
  the a, b and n for each processor. At that point it starts the timer for 
  the elapsed time and doesn't stop until the trapezoid rule has completed. 
  The function then adds up what each processor has calculated and prints out 
  all the results. 

  Adapted from Pacheco program 3.4 and program 5.1
*/
int main(void) {
  // Variable declarations
  int my_rank, comm_sz, n, local_n, src; 
  double true_error;
  double local_start, local_finish, local_elapsed, elapsed;
  double a, b, h, total, local_total, local_a, local_b;

  // Constant variables
  double true_val = 4003.7209001513268265;
  double target = .5 * pow(10, -12);

  // MPI declarations
  MPI_Init(NULL, NULL); 
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

  Get_input(my_rank, comm_sz, &a, &b, &n);

  h = (b-a)/n;  
  local_n = n/comm_sz;

  local_a = a + my_rank*local_n*h;
  local_b = local_a + local_n*h;

  MPI_Barrier(MPI_COMM_WORLD);
  if(my_rank == 0)
    local_start = MPI_Wtime();

  local_total = Trap(local_a, local_b, local_n, h);

  if(my_rank == 0) {
    local_finish = MPI_Wtime();
    local_elapsed = local_finish - local_start;
  }

  MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  // Add calculated totals from each processor
  if(my_rank != 0) {
    MPI_Send(&local_total, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
  else {   
    total = local_total;
    for (src = 1; src < comm_sz; src++) {
      MPI_Recv(&local_total,1, MPI_DOUBLE, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      total += local_total;
    }
  }

  // Print calculated data
  if(my_rank == 0) {
    printf("Running on %d processors.\n", comm_sz);
    printf("Elapsed time = %e seconds\n", elapsed);
    printf("With n = %d trapezoids, our estimate\n", n);
    printf("of the integral from %f to %f = %.13e\n", a, b, total);
    
    if(total - true_val < 0)
      true_error = (total - true_val)*-1;
    else
      true_error = total - true_val;
    
    if(true_error/true_val > target)
      printf("absolute relative true error = %e is NOT less than criteria %e \n", 
        true_error/true_val, target);
    else 
      printf("absolute relative true error = %e is less than criteria %e \n", 
        true_error/true_val, target);
    //printf("\n%.10f\n", h);
  }

  MPI_Finalize();

  return 0;
}