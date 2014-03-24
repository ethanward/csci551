#include <stdio.h>
#include <math.h>

double target = .5 * pow(10, -12);
double true_val = 4003.7209001513268265;



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
//2001875
//964104
int min_trap_count(double left, double right) {
  int count = 964100, min_count;
  double estimate, min = 1000000000000;

  while(count > 0) {
  	estimate = Trap(left, right, count, (right-left)/count);
    printf("%d:  %.14f \n", count, estimate);
    /*if(estimate < min) {
      min = estimate;
      min_count = count;
      printf("Min: %.14f\n", (min - true_val)/true_val);
      //printf("%d:  %.14f \n", count, (estimate - true_val)/true_val);
    }*/
  	if((estimate - true_val)/true_val < target) {
      printf("Estimate:  %.14e\n", estimate);
      printf("Error: %e\nCriteria: %e \n", (estimate - true_val)/true_val, target);
      return count;
    }
   
      //printf("%d:  %.10f \n", count, abs(estimate));
      count ++;
  }

  return -1;
}

int main() {
	int min = min_trap_count(100, 600);
	printf("Min: %d\n", min);
  //cout << "Min is: " << min << endl;
}
