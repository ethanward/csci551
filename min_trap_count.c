#include <stdio.h>
#include <math.h>


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
//2001875
int min_trap_count(double left, double right) {
  int count = 2001870, max_count = 1000000;
  double estimate;

  while(count > 0) {
  	estimate = Trap(left, right, count, (right-left)/count);
    printf("%d:  %.16f \n", count, estimate);
  	if(estimate - true_val < 0.0000000001) {
     // printf("%d:  %.11f\n", count, abs(true_val - estimate));
      printf("FUCK\n");
      return count;
    }
   
      //printf("%d:  %.10f \n", count, abs(estimate));
      count ++;
  }

  return -1;
}

int main() {
	int min = min_trap_count(100, 600);
	//cout << "Min is: " << min << endl;
}
