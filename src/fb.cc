#include "fb.h"

#include "matrix.h"

extern "C" {

// 	# Forward-Backward algorithm (used in Baum-Welch)
// 	# Returns alpha, beta, and full data likelihood
// 	
// 	# NOTE THE CHANGE IN FROM ROW TO COLUMN SUCH THAT TRANSPOSING A IS NOT NECCESSARY ANYMORE
// 	# IN COMPUTING ALPHA AND BETA BUT IS NOW NECCESSARY IN COMPUTING XI
// 	# A = T*K*K matrix with transition probabilities, from row to column!!!!!!!
// 	# B = T*K matrix with elements ab_{ij} = P(y_i|s_j)
// 	# init = K vector with initial probabilities
// 
// 	# NOTE: to prevent underflow, alpha and beta are scaled, using sca
// 	
// 	# NOTE: xi[t,i,j] = P(S[t] = j & S[t+1] = i) !!!NOTE the order of i and j!!!
	
// 	fb <- function(init,A,B,ntimes=NULL,return.all=FALSE,stationary=TRUE) {
// 
// 	fb(init=object@init,A=object@trDens,B=object@dens,ntimes=ntimes(object), 
// 			stationary=object@stationary,return.all=return.all)
	

// inputs are:
// a) ns: the number of states
// b) nc: the number of cases
// c) nt: the number of rows of data
// d) ntimes: rows of data of individual cases
// 1) init: ns vector or ns by nc matrix
// 2) trdens: ns by ns matrix or nt by ns by ns array
// 3) dens: nt by ns matrix
// 
// outputs are:
// 1) alpha: nt by ns matrix
// 2) beta: nt by ns matrix
// 3) xi: nt by ns by ns array
// 4) sca: nt vector

// gamma is computed as alpha*beta/sca in R (no loop needed)


void forwardbackward(int *ns, int *nc, int *nt, int *ntimes, 
					 double *init, double *trdens, double *dens, 
					 double *alpha, double *beta, double *sca, double *xi) {
		
		Rprintf("ns=%d\n",ns[0]);
		Rprintf("nc=%d\n",nc[0]);
		Rprintf("nt=%d\n",nt[0]);
		
		int nttot=0;
		
 		for(int i=0; i<3; i++) {
 			nttot += ntimes[i];
 			Rprintf("ntimes: %d \n", ntimes[i]);
 			Rprintf("nttot: %d \n", nttot);
 		}
		
		// 	loop over cases
		int tt=0;
		for(int cas=0; cas<nc[0]; cas++) {
			// compute alpha1 for this case
			double sca1=0.0;
			for(int i=0; i<ns[0]; i++) {
				alpha[tt*ns[0]+i] = init[cas*ns[0]+i]*dens[tt*ns[0]+i];
				sca1 += alpha[tt*ns[0]+i];
			}
			sca[tt] = sca1;
			for(int i=0; i<ns[0]; i++) {
				alpha[tt*ns[0]+i] /= sca1;
			}
			
			// compute scale for alpha1
			
			// scale alpha1
			
			//loop over ntimes[cas]>1
			
			// compute cumulative t for indexing alpha
			if(cas>0) tt += ntimes[cas];
			
			Rprintf("%d\n",tt);
			
		} // end cases
		
		
		// (auxiliaries for) gradients from mixing proportions
// 		matrix *xiC; xiC = new matrix[2];
// 		mp1[0](1)=3.7;
// 		mp1[1](1)=1.6;	
// 		mp1[0].print();
		
	}

	
// //this is the actual loglikelihood function for models with covariates
// void forwardbackward(double logl, double *init, int *linit) {
// 	
// 	double crit=0.0;
// 	
// 	crit=init[0];
// 	// This will stop optimization when the loglikelihood is better then crit, useful for
// 	// goodness of fit bootstrapping
// 	Rprintf("stop crit=%f\n",crit);
// 	
// 	Rprintf("Starting to compute likelihood.\n");
// 	
// 	Rprintf("linit: %d \n", linit[0]);
// 	
// 	logl=0;
// 	
// 	Rprintf("logl: %f \n", logl);
// 
// 	double tmp;
// 	
// 	for(int i=0; i<3; i++) {
// 		logl += init[i];
// 		tmp = init[i];
// 		Rprintf("init: %f \n", tmp);
// 		Rprintf("logl: %f \n", logl);
// 	}
// 	
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
		
// }
 


} // end extern "C"