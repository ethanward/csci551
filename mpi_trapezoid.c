#include <stdio.h>
#include "mpi.h"
#include <math.h>
//#include "trap.h"

// Approc Val: 4003.7209001513997464
// True Val: 4003.7209001513268265
// Abs Rel True Err: 1.82130327808923642875 × 10^-14

double true_val = 4003.7209001513;

double f(double x) { 
  return (cos(x/3) - 2*cos(x/5) + 5*sin(x/4) + 8); 
}

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
   } else { /* my_rank != 0 */
      MPI_Recv(a_p, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);
      MPI_Recv(b_p, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);
      MPI_Recv(n_p, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);
   } 
} 

int main(void) {
  int my_rank, comm_sz, n, local_n, src; 
  double true_val = 4003.7209001513268265;
  double target = .5 * pow(10, -12);
  double local_start, local_finish, local_elapsed, elapsed;
  double a, // The user input for a 
    local_total,
    total,
    b, // The user input for b
    h, // The user input for h
    local_a, // The calculated a for each proc
    local_b; // The calculated b for each proc 

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

  if(my_rank == 0) {
    printf("Running on %d processors.\n", comm_sz);
    printf("Elapsed time = %e seconds\n", elapsed);
    printf("With n = %d trapezoids, our estimate\n", n);
    printf("of the integral from %f to %f = %.13e\n", a, b, total);
    printf("absolute relative true error = %e is NOT less than criteria %e \n", (total - true_val)/true_val, target);
    //printf("\n%.10f\n", h);
  }

  MPI_Finalize();

  return 0;
}