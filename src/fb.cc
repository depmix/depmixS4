#include "fb.h"

extern "C" {
	
//this is the actual loglikelihood function for models with covariates
double forwardbackward(double *init, int *linit) {
	
	Rprintf("Starting to compute likelihood.\n");
	
	int ltin;
	
	ltin = linit[0];
	 
// 	Rprintf(ltin);
	
	double res;
	res=0;
	
	for(int i=0; i<ltin; i++) {
		res += init[i];
	}
	
	
// 	int nrcomp=models.getNrComp();
// 	int ngroups=models.getNrGroups();
// 	int *npars = new int[nrcomp];
// 	int *nstates = new int[nrcomp];
// 	int totalstates=0;
// 	
// 	for(int cp=0; cp<nrcomp; cp++) {
// 		npars[cp]=models.mods[cp].getNPars();
// 		nstates[cp]=models.mods[cp].getStates();
// 		totalstates += nstates[cp];
// 	}
// 	
// 	// (auxiliaries for) gradients from mixing proportions
// 	matrix *mp1; mp1 = new matrix[2];
// 	matrix *mpt; mpt = new matrix[ngroups];
// 	matrix *mptfinal; mptfinal = new matrix[ngroups];
	
	return(res);
	
} 


} // end extern "C"