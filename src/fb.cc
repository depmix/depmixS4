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


void forwardbackward(int *ns, int *nc, int *nt, int *ntimes, int *bt, int *et,
					 double *init, double *trdens, double *dens, 
					 double *alpha, double *beta, double *sca, double *xi) {
		
	// compute forward variables
	// loop over cases
	for(int cas=0; cas<nc[0]; cas++) {
		// compute alpha1 for this case
		double sca1=0.0;
		matrix alpha1(ns[0],1);
		matrix denst(ns[0]);
		for(int i=0; i<ns[0]; i++) {
			alpha1(i+1) = init[cas*ns[0]+i]*dens[(bt[cas]-1)*ns[0]+i];
			denst(i+1) = dens[(bt[cas]-1)*ns[0]+i];
			// compute scale for alpha1
			sca1 += alpha1(i+1);
		}
		sca[(bt[cas]-1)] = 1/sca1;
		// scale alpha1 and copy to alpha
		for(int i=0; i<ns[0]; i++) {
			alpha1(i+1) *= sca[(bt[cas]-1)];
			alpha[(bt[cas]-1)*ns[0]+i] = alpha1(i+1);
		}
		matrix alphat(ns[0],1);
		matrix trans(ns[0],ns[0]);
		// loop over ntimes[cas]>1
		if(ntimes[cas]>0) {
			for(int t=bt[cas]; t<et[cas]; t++) {
				// get trans and dens values
				for(int i=0; i<ns[0]; i++) {
					for(int j=0; j<ns[0]; j++) {
						trans(i+1,j+1) = trdens[(i*ns[0]+j)*nt[0]+t-1]; // A_t
					}
					denst(i+1) = dens[t*ns[0]+i]; //B_t+1
				}					
				// compute alphat
				alphat = had(transpose(trans)*alpha1,denst);// A_t*alpha_t*B_t+1
				// compute scale for t
				sca[t] = 1/(alphat.msum());
				// scale alphat values
				for(int i=0; i<ns[0]; i++) {
					alphat(i+1) *= sca[t];
					alpha[t*ns[0]+i] = alphat(i+1);
				}
				alpha1 = alphat;
			}
		}			
	} // end cases
	
	// compute backward variables and xi
	
	
}

} // end extern "C"