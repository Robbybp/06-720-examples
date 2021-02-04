#define MAX_ITER 50
#define TOLERANCE 1e-8
#define VERBOSE 0
#define DEBUG 0

/* ASL includes */
#include "asl.h"
#include "getstub.h"

#include <string.h>
#include <assert.h>

/* Local includes */
#include "convergence.h"
#include "linalg.h"
#include "iter_log.h"

char solver_name[] = "Newton solver";
static Option_Info Oinfo;
/* Ideally, we communicate solver options via an environment variable
 * or the command line. These options would be delivered to our solver
 * by Oinfo. We are not doing that here, but ASL still seems to require
 * this variable to get the nl file.
 */

int main(int argc, char **argv){
#if VERBOSE
  // First argument appears reserved for the program itself
  printf("Received %i argument(s)\n", argc-1);
  for (int i=1; i<argc; i++){
    printf("    Argument %i: %s\n", i, argv[i]);
  }
  printf("\n");
#endif

  ASL *asl;
  asl = ASL_alloc(ASL_read_fg);

  char sol_msg[100];

  char *stub = getstops(argv, &Oinfo);
  if (!stub) return 1;
  /* Not sure if getstops is necessary (as opposed to just getstub),
   * but want to follow the ASL docs as much as possible.
   * Note that argv here is type char**. This is consistent with
   * solvers/getstub.h.
   */
  printf("------------------------------------------\n");
  printf("Demo Newton solver. Written for 06-720 S21\n");
  printf("------------------------------------------\n\n");

  printf("Solving problem %s.nl\n", stub);

  FILE *fp = jac0dim(stub, (fint)strlen(stub));
  /* See ASL documentation.
   * Initializes basic info. For our purposes, after this call
   * the n_var and n_con variables have been initialized.
   */
  printf("n_var = %i\nn_con = %i\n", n_var, n_con);
  assert(n_var == n_con);

  double *x = (double *)malloc(sizeof(double) * n_var);
  double *x_plus = (double *)malloc(sizeof(double) * n_var);

  /* To store the value of the function f(x) = 0 */
  double *f = (double *)malloc(sizeof(double) * n_con);

  /* This tells asl to allocate initial primal array asl->i.X0_ */
  asl->i.want_xpi0_ = 1;
  /* This reads expressions from the nl file */
  fg_read(fp, 0);

#if VERBOSE
  printf("\nAt the initial guess:\n");
#endif
  /* Initialize our local "current x" array */
  for (int i=0; i<n_var; i++){
    x[i] = asl->i.X0_[i];
    x_plus[i] = x[i];
#if VERBOSE
    printf("x_%i: %f\n", i, x[i]);
#endif
  }

  conval(x, f, 0);
  assert(!Urhsx);
  /* ^ This is some behavior of ASL I don't really understand.
   * It appears the bounds can be stored either in the same array
   * or in separate arrays. The default appears to be the former but
   * I don't see how or why you might change this.
   */
  for (int j=0; j<n_con; j++){
    /* Require all constraints to be equality constraints */
    assert(LUrhs[2*j] == LUrhs[2*j+1]);
#if VERBOSE
    printf("%f <= (f_%i = %f) <= %f\n", LUrhs[2*j], j, f[j], LUrhs[2*j+1]);
#endif
    f[j] -= LUrhs[2*j];
  }

  // asl.h defines: nzc = asl->i.nzc_
  double *J = (double *)malloc(sizeof(double)*nzc);

  /* Allocate arrays to store COO repr of Jacobian */
  int *row = (int *)malloc(sizeof(double)*nzc);
  int *col = (int *)malloc(sizeof(double)*nzc);
  double *entries = (double *)malloc(sizeof(double)*nzc);

  /* Search direction */
  double *d = (double *)malloc(sizeof(double)*n_var);

  /* -------------- */
  /* Algorithm loop */
  /* -------------- */
  int iter = 0;
  print_header();
  while (iter < MAX_ITER){
    double eps_conv = 0.5*inner_product(n_con, f, f);
    log_iter(iter, eps_conv);
    if (check_conv_ip(n_con, f, TOLERANCE)) {
      printf("\nSolution found.\n||f|| = %e\n", eps_conv);
#if VERBOSE
      printf("\nVariable values are:\n");
      for (int i=0; i<n_var; i++){
        printf("x_%i: %f\n", i, x[i]);
      }
#endif
      strcpy(sol_msg, "Solution found.");
      solve_result_num = 0;
      Oinfo.wantsol = 9;
      /* ^ (write a sol file & suppress printing of solution message) */
      write_sol(sol_msg, x_plus, NULL, &Oinfo);
      break;
    }

    // Advance x after checking convergece
    for (int i=0; i<n_var; i++){
      x[i] = x_plus[i];
    }

    /* Note that conval has already been called with the appropriate
     * variable vector.
     */
    jacval(x, J, 0);

    /* Now populate COO arrays from ASL's linked-list data structure.
     * (The data structure somewhat approximates CSR)
     */
    cgrad *cg;
    int k = 0;
    for (int j=0; j<n_con; j++){
      cg = Cgrad[j];
      do{
        // Assuming that every constraint has a variable here...
        row[k] = j;
        col[k] = cg->varno;
        entries[k] = J[cg->goff];
        k++;
      }
      while((cg = cg->next));
    }
    assert(k == nzc);

#if DEBUG
    for (int k=0; k<nzc; k++){
      printf("%i, %i, %f\n", row[k], col[k], entries[k]);
    }
#endif

    int status = solve_with_ma28(n_var, nzc, row, col, entries, f, d);
#if DEBUG
    printf("MA28 exit status: %i\n", status);
#endif
    if (status < 0) {
      printf("Error in MA28. Aborting.\n");
      return -1;
    }
  
    /* We have solved J^-1 f = d.
     * Reverse sign of d to get Newton step -J^-1 f.
     */
    for (int i=0; i<n_var; i++){
      d[i] = -d[i];
    }

#if DEBUG
    for (int i=0; i<n_var; i++){
      printf("d_%i = %f\n", i, d[i]);
    }
#endif
  
    /* Update variable vector */
    for (int i=0; i<n_var; i++){
      x_plus[i] = x[i] + d[i];
    }

    /* Update f. Want to do this before the next time we check
     * convergence. In the future our globalization strategy may
     * depend on this value.
     */
    conval(x_plus, f, 0);
    /* Need to adjust f s.t. constraint is f(x) == 0
     */
    for (int j=0; j<n_con; j++){
      f[j] -= LUrhs[2*j];
    }

    iter += 1;

    //break;
  } // while(iter<MAX_ITER)


  /* TODO:
   * Add failure message if iter == MAX_ITER
   * Write sol file.
   * Function for printing useful information at every iteration.
   *
   * See compile.sh
   */

  free(entries);
  free(row);
  free(col);
  free(d);
  free(x);
  free(x_plus);
  free(f);
  ASL_free(&asl); // Takes type (ASL **) for some reason...
  return 0;
}
