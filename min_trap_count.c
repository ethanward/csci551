/*

  Author: Ethan Ward
  ID: 005704518

  This program is for finding the minimum number of trapezoids that 
  gives a correct answer in 14 significant digits for the definite
  integral of cos(x/3) - 2*cos(x/5) + 5*sin(x/4) + 8

*/

#include <stdio.h>
#include <math.h>

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
  Input: a and b for the trapezoidal rule
  Output: the minimum number of trapezoids that 
    gives a correct answer in 14 significant digits for the definite
    integral of cos(x/3) - 2*cos(x/5) + 5*sin(x/4) + 8
*/
int min_trap_count(double left, double right) {
  // Constant vaiables
  double target = 0.5 * pow(10, -12);
  double true_val = 4003.7209001513268265;
  int count = 964100; // 964104 is the true answer. This was used to
                      // check correctness

  // Variable declaration
  double estimate;

  // Essentially an infinite loop. Stops when absolute relative true error
  // is less than target 
  while(count > 0) {
    estimate = Trap(left, right, count, (right-left)/count);
    printf("%d:  %.14f \n", count, estimate);
    
    if((estimate - true_val)/true_val < target) {
      printf("Estimate:  %.14e\n", estimate);
      printf("Error: %e\nCriteria: %e \n", (estimate - true_val)/true_val, target);
      return count;
    }
    else
      count ++;
  }

  return -1;
}

int main() {
	int min = min_trap_count(100, 600);
	printf("Min: %d\n", min);

}
