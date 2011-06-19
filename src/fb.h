
#ifndef FB
#define FB 1

#include <stdio.h>
#include <stdlib.h>
/* #include <fstream.h> */

 
/* #include "matrix.h" */
  
extern "C" {
	
#include <R.h>    
#include <Rmath.h>

// criterion for stopping optimization, used in bootstrapping
void forwardbackward(int *ns, int *nc, int *nt, int *ntimes, double *init, double *trdens, double *dens, double *alpha, double *beta, double *sca, double *xi);

// the main function that computes the forward backward vars, and gamma and xi
/* void forwardbackward(double logl, double *init, int *linit); */
	
} //end extern "C"

#endif