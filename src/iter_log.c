#include <stdio.h>
#include "iter_log.h"

int print_header(void){
  printf("\n\tIter\t||f||\n");
  return 0;
}

int log_iter(int iter, double eps){
  printf("\t%i\t%e\n", iter, eps);
  return 0;
}
