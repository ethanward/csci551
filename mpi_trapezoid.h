#ifndef _mpi_trapezoid_h_
#define _mpi_trapezoid_h_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

class
	public:
		double Trap(double left, double right, int trap_count, double base);
		double f(double x);
		void Get_input(int my_rank, int comm_sz, double* a_p, double* b_p, int* n_p);
#endif