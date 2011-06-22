#include "ddens.h"

extern "C" {

void ddens(double *dens, double *res, int *dim) {
	
// 	Rprintf("dim 0: %d \n", dim[0]);
// 	Rprintf("dim 1: %d \n", dim[1]);
// 	Rprintf("dim 2: %d \n", dim[2]);
	
	for(int t=0; t<dim[0]; t++) {
// 		Rprintf("t: %d \n", t);
		for(int i=0; i<dim[1]; i++) {
// 			Rprintf("i: %d \n", i);
			res[t*dim[2]+i] = dens[i*dim[0]*dim[2]+t];
// 			Rprintf("dens %f \n", dens[i*dim[0]*dim[2]+t]);
			for(int j=1; j<dim[2]; j++) {
// 				Rprintf("j: %d \n", j);
				res[t*dim[2]+i] *= dens[i*dim[0]*dim[2]+j*dim[0]+t];
// 				Rprintf("dens %f \n", dens[i*dim[0]*dim[2]+j*dim[0]+t]);
				
			}
		}
	}	
}

} // end extern "C"