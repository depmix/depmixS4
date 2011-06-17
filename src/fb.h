
#include <stdio.h>
#include <stdlib.h>
/* #include <fstream.h> */

#include "matrix.h"

extern "C" {

#include <R.h>	
#include <Rmath.h>

// the main function that computes the forward backward vars, and gamma and xi
double forwardbackward(double *init, int *linit);

} //end extern "C"

