#include "convergence.h"

/* Implementations of functions used to check convergence.
 */

/* Computes the standard inner product between two vectors of
 * the same dimension.
 */
double inner_product(int dim, double *v1, double *v2){
  double ip = 0.0;
  for (int i=0; i<dim; i++){
    ip += v1[i]*v2[i];
  }
  return ip;
}

/* Checks whether a vector is zero by comparing one half the
 * inner product to some tolerance.
 */
int check_conv_ip(int dim, double *f, double tol){
  double ip = inner_product(dim, f, f);
  return (int)(0.5*ip <= tol);
}
