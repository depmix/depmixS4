em <- function(object,maxit=100,tol=1e-6,verbose=FALSE,...) {
	
	if(!is(object,"depmix")) stop("object is not of class 'depmix'")
	
	ns <- object@nstates
	
	ntimes <- ntimes(object)
	lt <- length(ntimes)
	et <- cumsum(ntimes)
	bt <- c(1,et[-lt]+1)
		
	LL <- logLik(object)
	
	converge <- FALSE
	j <- 0
	
	A <- object@trDens
	while(j <= maxit & !converge) {
		
		for(i in 1:ns) {
			A[,,i] <- object@trDens[,,i]
		}
		
		B <- apply(object@dens,c(1,3),prod)
		init <- object@init
		LL.old <- LL
		j <- j+1
		
		# expectation
		fbo <- fb(init=init,A=A,B=B,ntimes=ntimes(object))
		LL <- fbo$logLike
		
		# maximization
				
		# should become object@prior <- fit(object@prior)
		object@prior@y <- fbo$gamma[bt,,drop=FALSE]
		object@prior <- fit(object@prior, w=NULL,ntimes=NULL)
		object@init <- dens(object@prior)
				
		trm <- matrix(0,ns,ns)
		for(i in 1:ns) {
			
			if(max(ntimes(object)>1)) { # skip transition parameters update in case of latent class model
				if(!object@stationary) {
					object@transition[[i]]@y <- fbo$xi[,,i]/fbo$gamma[,i]
					object@transition[[i]] <- fit(object@transition[[i]],w=as.matrix(fbo$gamma[,i]),ntimes=ntimes(object)) # check this
				} else {
					for(k in 1:ns) {
						trm[i,k] <- sum(fbo$xi[-c(et),k,i])/sum(fbo$gamma[-c(et),i])
					}
					# FIX THIS; it will only work with a specific trinModel
					object@transition[[i]]@parameters$coefficients <- object@transition[[i]]@family$linkfun(trm[i,],base=object@transition[[i]]@family$base)
				}
				
				# update trDens slot of the model
				object@trDens[,,i] <- dens(object@transition[[i]])
			}
			
			for(k in 1:nresp(object)) {
				object@response[[i]][[k]] <- fit(object@response[[i]][[k]],w=fbo$gamma[,i])
				# update dens slot of the model
				object@dens[,k,i] <- dens(object@response[[i]][[k]])
			}
		}
				
		LL <- logLik(object)
		if(verbose&((j%%5)==0)) cat("iteration",j,"logLik:",LL,"\n")
		if( (LL >= LL.old) & (LL - LL.old < tol))  {
			cat("iteration",j,"logLik:",LL,"\n")
			converge <- TRUE
		}
	}
	
	class(object) <- "depmix.fitted"
	if(converge) object@message <- "Log likelihood converged to within tol."
	else object@message <- "'maxit' iterations reached in EM without convergence."
	
	# no constraints in EM
	object@conMat <- matrix()
	
	# what do we want in slot posterior?
	object@posterior <- fbo$gamma
	
	object
}
