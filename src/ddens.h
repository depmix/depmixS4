
#ifndef DDENS
#define DDENS 1

#include <stdio.h>
#include <stdlib.h>
  
extern "C" {
	
#include <R.h>    
#include <Rmath.h>

// compute the product of densities for each t
// hopefully more efficient than apply?!?!?
void ddens(double *dens, double *densOut, int *dim);
	
} //end extern "C"

#endif