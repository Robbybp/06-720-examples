#include <stdio.h>
#include "linalg.h"
#include "asl.h"

int solve_with_ma28(fint dim, fint nnz, fint *row, fint *col,
    double *A, double *b, double *x){
  // Not immediately clear when I should use hsl_int vs fint

  double U = 0.1;
  hsl_int *IKEEP = (hsl_int *)malloc(sizeof(hsl_int)*5*dim);
  hsl_int *IW = (hsl_int *)malloc(sizeof(hsl_int)*8*dim);
  double *W = (double *)malloc(sizeof(double)*dim);
  int IFLAG;
  int MTYPE = 1;

  /* Populate x with entries from RHS b. These will be overridden
   * by the call to ma28c.
   */
  for (int i=0; i<dim; i++){
    x[i] = b[i];
  }

  /* Convert row and col arrays to be compatible with Fortran and 
   * allocate enough memory to store the factors.
   * Estimates for necessary storage space are taken from the
   * documentation.
   */
  int LICN = 4*nnz;
  int LIRN = 2*nnz;
  int *IRN = (int *)malloc(sizeof(int)*LIRN);
  int *JCN = (int *)malloc(sizeof(int)*LICN);
  double *A_factor = (double *)malloc(sizeof(double)*LICN);
  for (int k=0; k<nnz; k++){
    IRN[k] = row[k] + 1;
    JCN[k] = col[k] + 1;
    A_factor[k] = A[k];
  }

  ma28ad_(&dim, &nnz, A_factor, &LICN, IRN, &LIRN, JCN, &U, IKEEP, IW, W, &IFLAG);
  ma28cd_(&dim, A_factor, &LICN, JCN, IKEEP, x, W, &MTYPE);

  free(A_factor);
  free(IRN);
  free(JCN);
  free(W);
  free(IW);
  free(IKEEP);
  return IFLAG;
}
