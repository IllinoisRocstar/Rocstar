/*
// serialBS95.c
//
// Functions to be called by FORTRAN that will solve sparse 
// linear systems of equations using BlockSolve95.
//
//
// NOTE:  FLOAT is defined as double, not float
//
//
// TO RUN FOLLOW THESE STEPS:
//  1 - Set up the problem
//  2 - Run BS95SETUP
//  3 - Run BS95SOLVE
//  4 - Repeat step 4 as needed if system needs to be solved with different rhs vectors
//  5 - Run BS95FREE
//  6 - Repeat steps 1-5 if system needs to be solved with different A matrices
//  7 - Run BS95FINALIZE to finalize MPI communications
//
*/


/*
// Include standard and BlockSolve95 headers
*/
#include <stdlib.h>
#include <stdio.h>

#include <BSprivate.h>

//#include "mpi.h"
extern MPI_Comm GetCommunicator();
/*
// Global variables
*/
BSprocinfo *procinfo;
BSspmat *bA;
BSpar_mat *pA, *f_pA;
BScomm *Acomm, *fcomm;
int maxiter = 10000;




/*
// Utility functions
*/



/*
 * NAME
 *    ctranspose
 *
 * FUNCTION
 *    Transposes a matrix from FORTRAN style
 *    to C style.
 *
 * INPUTS
 *    n -- size of one dimension of the matrix
 *    d -- the square matrix to be transposed
 *
 * OUTPUTS
 *    d[n][n]  -- the transposed square matrix
 *
 * USES
 *    none
 *
 */

void ctranspose(int n,double d[n][n]) {
  double tempdmat[n][n];
  int i;
  int j;
  for(i=0;i<n;i++) {
     for(j=0;j<n;j++) {
       tempdmat[i][j] = d[j][i];
     }
  }
  for(i=0;i<n;i++) {
     for(j=0;j<n;j++) {
       d[i][j] = tempdmat[i][j];
     }
   }
}




/*
// Functions to be called by Fortran
*/





/*
 * NAME
 *    bs95setup
 *
 * FUNCTION
 *    Sets up everything needed to solve a system of
 *    linear equations in parallel using BlockSolve95.
 *
 * INPUTS
 *    *ndim        -- dimension of the A matrix were it to be explicitly formed
 *    *nnz         -- number of non-zeros in the A matrix
 *    *nstart      -- first row that is assigned to this processor
 *    *nrows       -- number of rows assigned to this processor
 *    rp           -- mapping of positions in cval vector to rows (see BS95 manual for details)
 *    cval         -- mapping of values in aval to collums (see BS95 manual for details)
 *    aval         -- values of all the nonzeros in the A matrix (see BS95 manual for details)
 *    *symm        -- 1 if matrix nonzeros are symmetric, 0 otherwise
 *    *debugwrite  -- 1 if debug statements should be written to the screen, 0 otherwise
 *
 * OUTPUTS
 *    none
 *
 * USES
 *    BlockSolve95
 *
 */

void bs95setup_(int *ndim, int *nnz, int *nstart, int *nrows, int rp[*nrows+1], int cval[*nnz], double aval[*nnz], int *symm, int *debugwrite) {

  int argc = 1;
  int i;
  char **argv;
  float my_alpha;

  if(*debugwrite==1) {
    printf("ndim = %d\n",*ndim);
    printf("nnz = %d\n",*nnz);
    printf("nstart = %d\n",*nstart);
    printf("nrows = %d\n",*nrows);
    printf("rp = [ ");
    for(i=0;i<*nrows+1;i++) {
      printf("%d ",rp[i]);
    }
    printf("]\ncval =  [ ");
    for(i=0;i<*nnz;i++) {
      printf("%d ",cval[i]);
    }
    printf("]\naval =  [ ");
    for(i=0;i<*nnz;i++) {
      printf("%f ",aval[i]);
    }
    printf("]\n");
  }
  
  
  /*
  // Initialize BlockSolve with dummy command-line arguments
  */
  argv = (char ** ) MALLOC(sizeof(char)*argc);
  for(i=0;i<argc;i++) { 
    argv[i] = "1";
  }
  if(*debugwrite==1){printf("   BS95SETUP\n     initializing BS...");}
  BSinit(&argc,&argv);
  if(*debugwrite==1){printf("done\n");}

  /*
  // Get the processor context
  */
  if(*debugwrite==1){printf("     getting the processor context...");}
  procinfo = BScreate_ctx();
  //  MPI_Comm self = MPI_COMM_SELF;
  MPI_Comm the_communicator = GetCommunicator();
  BSctx_set_ps(procinfo,the_communicator);
  BSctx_set_err(procinfo,1);
  if(*debugwrite==1){printf("done\n");}

  /*
  // Create the A matrix in BSspmat datatype
  */
  if(*debugwrite==1){printf("     creating the BSspmat A matrix...");}
  bA = BSeasy_A(*nstart,*nrows,rp,cval,aval,procinfo);
  if(*debugwrite==1){printf("done\n");}

  /*
  // Define whether matrix nonzeros are symmetric
  */
  if(*debugwrite==1){printf("     setting up matrix storage...");}
  BSset_mat_symmetric(bA,*symm);
  if(*debugwrite==1){printf("done\n");}

  /*
  // Set up the storage type (Cholesky or ILU)
  */
  if(*debugwrite==1){printf("     setting up storage type...");}
  BSset_mat_icc_storage(bA,*symm);
  if(*debugwrite==1){printf("done\n");}

  /*
  // Convert matrix to BSpar_mat datatype
  */
  if(*debugwrite==1){printf("     permuting the A matrix...");}
  pA = BSmain_perm(procinfo,bA);
  if(*debugwrite==1){printf("done\n");}

  /*
  // Diagonally scale the matrix
  */
  if(*debugwrite==1){printf("     diagonally scaling the A matrix...");}
  BSscale_diag(pA,pA->diag,procinfo);
  if(*debugwrite==1){printf("done\n");}

  /*
  // Set up communication data structures
  */
  if(*debugwrite==1){printf("     setting up communication data structures...");}
  Acomm = BSsetup_forward(pA,procinfo);
  fcomm = BSsetup_factor(pA,procinfo);
  if(*debugwrite==1){printf("done\n");}

  /*
  // Copy the pA matrix
  */
  if(*debugwrite==1){printf("     creating a copy of the permuted matrix...");}
  f_pA = BScopy_par_mat(pA);
  if(*debugwrite==1){printf("done\n");}

  /*
  // Factor the matrix
  */
  if(*debugwrite==1){printf("     performing an incomplete factorization of the matrix...");}
  my_alpha = 1.0;
  while (BSfactor(f_pA,fcomm,procinfo) != 0) {
    BScopy_nz(pA,f_pA);
    my_alpha += 0.1;
    BSset_diag(f_pA,my_alpha,procinfo);
  }
  if(*debugwrite==1){printf("done\n");}

  /*
  // Set solver parameters
  */
  if(*debugwrite==1){printf("     setting solver parameters...");}
  if(*symm==1) { 
    /* // Use CG if the system is symmetric */
    BSctx_set_pre(procinfo,PRE_ICC);
    BSctx_set_method(procinfo,CG);
  }
  else { 
    /* // Use GMRES if the system is non-symmetric */
    BSctx_set_pre(procinfo,PRE_ILU);
    BSctx_set_method(procinfo,GMRES);
  }
  BSctx_set_max_it(procinfo,maxiter);
  BSctx_set_guess(procinfo,FALSE);  /* // Use vector passed in as initial guess */
  if(*debugwrite==1){printf("done\n");}

}








/*
 * NAME
 *    bs95solve
 *
 * FUNCTION
 *    Solves a set of linear equations using BlockSolve95.
 *    Previously, bs95setup must have been run in order to
 *    set up the A matrix.
 *
 * INPUTS
 *    *ndim        -- dimension of the rhs and x vectors
 *    rhs          -- the right-hand side of the A*x=rhs equation
 *    x            -- the solution to the A*x=rhs equation
 *    *debugwrite  -- 1 if debug statements should be written to the screen, 0 otherwise
 *
 * OUTPUTS
 *    none
 *
 * USES
 *    BlockSolve95
 *
 */


void bs95solve_(int *ndim, double rhs[*ndim], double x[*ndim], double *tolerance, int *debugwrite) {

  FLOAT residual;
  int num_iter;
  int i;

  if(*debugwrite==1) {
    printf("rhs = [ ");
    for(i=0;i<*ndim;i++) {
      printf("%f ",rhs[i]);
    }
    printf("]\n");
  }

  /*
  // Solve the system
  */
  if(*debugwrite==1){printf("   BS95SOLVE\n     solving the system...");}
  BSctx_set_tol(procinfo,*tolerance);
  num_iter = BSpar_solve(pA,f_pA,Acomm,rhs,x,&residual,procinfo);
  if(num_iter >= maxiter) {
    printf("*** WARNING ***  BlockSolve95 using too many iterations.  Matrix likely ill-conditioned.\n");
  }
  if(*debugwrite==1){
    printf("done\n");
    printf("       BS took %d iterations to solve the system\n",num_iter);
    printf("       x = [ ");
    for(i=0;i<*ndim;i++) { printf("%f ",x[i]); }
    printf("]\n");
    printf("       residual = %e\n",residual);
  }
  /* //printf("       BS took %d iterations to solve the system\n",num_iter); */

}





/*
 * NAME
 *    bs95free
 *
 * FUNCTION
 *    Frees the contexts and memory used by BlockSolve95 when
 *    it solves a set of linear equations.
 *
 * INPUTS
 *    *debugwrite  -- 1 if debug statements should be written to the screen, 0 otherwise
 *
 * OUTPUTS
 *    none
 *
 * USES
 *    BlockSolve95
 *
 */

void bs95free_(int *debugwrite) {

  /*
  // Free matrices
  */
  if(*debugwrite==1){printf("   BS95FREE\n     freeing matrices...");}
  BSfree_par_mat(pA);
  BSfree_copy_par_mat(f_pA);
  BSfree_easymat(bA);
  if(*debugwrite==1){printf("done\n");}

  /*
  // Free the communication patterns
  */
  if(*debugwrite==1){printf("     freeing the communication patterns...");}
  BSfree_comm(Acomm);
  BSfree_comm(fcomm);
  if(*debugwrite==1){printf("done\n");}

  /*
  // Free the processor context
  */
  if(*debugwrite==1){printf("     freeing the processor context...");}
  BSfree_ctx(procinfo);
  if(*debugwrite==1){printf("done\n");}

}




/*
 * NAME
 *    bs95finalize
 *
 * FUNCTION
 *    Finalizes BlockSolve95.  This should be only run once you
 *    are completely finished using BlockSolve95
 *
 * INPUTS
 *    *debugwrite  -- 1 if debug statements should be written to the screen, 0 otherwise
 *
 * OUTPUTS
 *    none
 *
 * USES
 *    BlockSolve95
 *
 */

void bs95finalize_(int *debugwrite) {

  /*
  // Finalize BlockSolve
  */
  if(*debugwrite==1){printf("   BS95FINALIZE\n     finalizing BS...");}
  BSfinalize();
  if(*debugwrite==1){printf("done\n");}

}




/*
// Generalized serial functions to be called by Fortran
*/




/*
 * NAME
 *    serialbs95tridiag
 *
 * FUNCTION
 *    Used to solve a tridiagonal matrix in serial (1 processor) using
 *    BlockSolve95 functions.
 *
 * INPUTS
 *    *ndim -- The dimension of one side of the square A matrix
 *    ad -- Vector containing lower diagonal entries of the matrix
 *    bd -- Vector containing on-diagonal entries of the matrix
 *    cd -- Vector containing upper diagonal entries of the matrix
 *    rhs -- The right-hand side vector
 *    *debugwrite -- 1 if debug statements should be written to the screen, 0 otherwise
 *
 * OUTPUTS
 *    x -- The solution vector
 *
 * USES
 *    BlockSolve95
 *
 */

void serialbs95tridiag_(int *ndim, double ad[*ndim], double bd[*ndim], double cd[*ndim], double rhs[*ndim], double x[*ndim], int *debugwrite) {

  /*
  // Solve tridiagonal systems of equations:
  //
  //     [ bd[0]     cd[0]      0         ...  0       ]   [   x[0]    ]   [    rhs[0]   ]
  //     [ ad[1]     bd[1]      cd[1]  0  ...  0       ]   [   x[1]    ]   [    rhs[1]   ]
  //     [      .         .          .                 ]   [     .     ]   [      .      ]
  //     [          .         .          .             ] * [     .     ] = [      .      ]
  //     [             .          .          .         ]   [     .     ]   [      .      ]
  //     [  0 ... 0  ad[ndim-2] bd[ndim-2] cd[ndim-2]  ]   [ x[ndim-2] ]   [ rhs[ndim-2] ]
  //     [  0 ...       0       ad[ndim-1] bd[ndim-1]  ]   [ x[ndim-1] ]   [ rhs[ndim-1] ]
  //
  */

  int *rp;
  int *cval;
  FLOAT *aval;
  int i;
  int nnz;
  int avali;
  int nstart = 0;
  int symm = 1;

  /*
  // Print header
  */
  if(*debugwrite==1){printf("   Using BlockSolve95 to solve the tridiagonal system...\n");}

  /*
  // Create the mapping of the nonzero entries of A
  */
  if(*debugwrite==1){printf("     defining matrix mapping...");}
  nnz = 3*(*ndim)-2;
  aval = (FLOAT *) MALLOC(sizeof(FLOAT)*nnz);
  cval = (int *) MALLOC(sizeof(int)*nnz);
  rp = (int *) MALLOC(sizeof(int)*(*ndim+1));
  printf(" \b");  /* // Give the program time to realize that new stuff is allocated */
  aval[0] = bd[0];
  aval[1] = cd[0];
  cval[0] = 0;
  cval[1] = 1;
  rp[0] = 0;
  rp[1] = 2;
  avali = 2;
  for(i=1;i<*ndim-1;i++) {
    aval[avali]   = ad[i];
    aval[avali+1] = bd[i];
    aval[avali+2] = cd[i];
    cval[avali]   = i-1;
    cval[avali+1] = i;
    cval[avali+2] = i+1;
    avali = avali + 3;
    rp[i+1] = avali;
  }
  aval[nnz-2] = ad[*ndim-1];
  aval[nnz-1] = bd[*ndim-1];
  cval[nnz-2] = *ndim-2;
  cval[nnz-1] = *ndim-1;
  rp[*ndim] = nnz;
  if(*debugwrite==1){printf("done\n");}
  
  /*
  // Solve the system
  */
  double *tol;
  *tol = 1.0e-09;
  bs95setup_(ndim,&nnz,&nstart,ndim,rp,cval,aval,&symm,debugwrite);
  bs95solve_(ndim,rhs,x,tol,debugwrite);
  bs95free_(debugwrite);

  /*
  // Free the matrix mapping
  */
  if(*debugwrite==1){printf("     freeing matrix mapping...");}
  free(aval);
  free(cval);
  free(rp);
  if(*debugwrite==1){printf("done\n");}

  /*
  // Print footer
  */
  if(*debugwrite==1){printf("   Done using BlockSolve95 to solve the system\n");}

}




/*
 * NAME
 *    serialbs95
 *
 * FUNCTION
 *    Used to solve a general system of linear equations serial  using
 *    BlockSolve95 functions.
 *
 * INPUTS
 *    *ndim -- The dimension of one side of the square A matrix
 *    A -- The matrix.
 *    rhs -- The right-hand side vector
 *    *debugwrite -- 1 if debug statements should be written to the screen, 0 otherwise
 *
 * OUTPUTS
 *    x -- The solution vector
 *
 * USES
 *    BlockSolve95
 *    ctranspose
 *
 */

void serialbs95_(int *ndim, double A[*ndim][*ndim], double rhs[*ndim], double x[*ndim], int *debugwrite) {

  /*
  // Solve general systems of equations
  //
  //     [A]*[x] = [rhs]
  //
  */

  int *rp;
  int *cval;
  FLOAT *aval;
  int i;
  int j;
  int nnz;
  int avali;
  int nstart = 0;
  int symm = 0;

  /*
  // Print footer
  */
  if(*debugwrite==1){printf("   Using BlockSolve95 to solve the system\n");}

  /*
  // Transpose the matrix indices of the 2D matrices from Fortran style to C style
  */
  ctranspose(*ndim,A);

  /*
  // Create the mapping of the nonzero entries of A
  */
  if(*debugwrite==1){printf("     defining matrix mapping...");}
  nnz = 0;
  for(i=0;i<*ndim;i++) {
    for(j=0;j<*ndim;j++) {
      if(A[i][j] != 0.0) {
	nnz++;
      }
    }
  }
  aval = (FLOAT *) MALLOC(sizeof(FLOAT)*nnz);
  cval = (int *) MALLOC(sizeof(int)*nnz);
  rp = (int *) MALLOC(sizeof(int)*(*ndim+1));
  printf(" \b");  /* // Give the program time to realize that new stuff has been allocated */
  avali = 0;
  rp[0] = 0;
  for(i=0;i<*ndim;i++) {
    for(j=0;j<*ndim;j++) {
      if(A[i][j] != 0.0) {
	aval[avali] = A[i][j];
	cval[avali] = j;
	avali++;
      }
    }
    rp[i+1] = avali;  }
  if(*debugwrite==1){printf("done\n");}

  /*
  // Solve the system
  */
  double *tol;
  *tol = 1.0e-09;
  bs95setup_(ndim,&nnz,&nstart,ndim,rp,cval,aval,&symm,debugwrite);
  bs95solve_(ndim,rhs,x,tol,debugwrite);
  bs95free_(debugwrite);

  /*
  // Free the matrix mapping
  */
  if(*debugwrite==1){printf("     freeing matrix mapping...");}
  free(aval);
  free(cval);
  free(rp);
  if(*debugwrite==1){printf("done\n");}

  /*
  // Transpose the matrix indices of the 2D matrices from C style to Fortran style
  */
  ctranspose(*ndim,A);

  /*
  // Print footer
  */
  if(*debugwrite==1){printf("   Done using BlockSolve95 to solve the system\n");}

}




/*
 * NAME
 *    serialbs95mat
 *
 * FUNCTION
 *    Takes a general matrix and puts it into compressed row storage (CRS)
 *    format for use with BlockSolve95 functions.
 *
 * INPUTS
 *    *ndim -- The dimension of one side of the square A matrix
 *    *nnz -- Reserved
 *    A -- The matrix.
 *    *debugwrite -- 1 if debug statements should be written to the screen, 0 otherwise
 *
 * OUTPUTS
 *    rp -- the row mapping CRS vector
 *    cval -- the collum mapping CRS vector
 *    aval -- the nonzero value CRS vector
 *
 * USES
 *    ctranspose
 *
 */

void serialbs95mat_(int *ndim, int *nnz, double A[*ndim][*ndim], int rp[*ndim+1], int cval[*nnz], double aval[*nnz], int *debugwrite) {

  /*
  // Construct the compressed row storage for a matrix
  */

  int avali;
  int i;
  int j;

  /*
  // Transpose the matrix indices of the 2D matrices from Fortran style to C style
  */
  ctranspose(*ndim,A);

  /*
  // Create the mapping of the nonzero entries of A
  */
  if(*debugwrite==1){printf("   SERIALBS95MAT\n     defining matrix mapping...");}
  avali = 0;
  rp[0] = 0;
  for(i=0;i<*ndim;i++) {
    for(j=0;j<*ndim;j++) {
      if(A[i][j] != 0.0) {
	aval[avali] = A[i][j];
	cval[avali] = j;
	avali++;
      }
    }
    rp[i+1] = avali;  
  }
  if(*debugwrite==1){printf("done\n");}

  /*
  // Transpose the matrix indices of the 2D matrices from C style to Fortran style
  */
  ctranspose(*ndim,A);

}
