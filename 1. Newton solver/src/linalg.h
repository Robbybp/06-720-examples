#include "asl.h"

/* ASL defines a "Fortran int" (long) that we're supposed to use to
 * as arguments for Fortran functions. Not sure if this is necessary,
 * so I typedef here so it can be changed easily for experimentation.
 */
typedef fint hsl_int;

/* External HSL functions */
/* Note that all variables are "pass by reference" in Fortran */
extern void ma28ad_(hsl_int *N, hsl_int *NZ, double *A, hsl_int *LICN,
    hsl_int *IRN, hsl_int *LIRN, hsl_int *ICN, double *U, hsl_int *IKEEP,
    hsl_int *IW, double *W, hsl_int *IFLAG);

extern void ma28bd_(hsl_int *N, hsl_int *NZ, double *A, hsl_int *LICN,
    hsl_int *IVECT, hsl_int *JVECT, hsl_int *ICN, hsl_int *IKEEP,
    hsl_int *IW, double *W, hsl_int *IFLAG);

extern void ma28cd_(hsl_int *N, double *A, hsl_int *LICN, hsl_int *ICN,
    hsl_int *IKEEP, double *RHS, double *W, hsl_int *MTYPE);

extern void ma28id_(hsl_int *N, hsl_int *NZ, double *AORG, hsl_int *IRNORG,
    hsl_int *ICNORG, hsl_int *LICN, double *A, hsl_int *ICN, hsl_int *IKEEP,
    double *RHS, double *X, double *R, double *W, hsl_int *MTYPE, double *PREC,
    hsl_int *IFLAG);

/* Functions to actually be used by solver */
int solve_with_ma28(fint dim, fint nnz, fint *row, fint *col,
    double *A, double *b, double *x);
