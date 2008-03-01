em <- function(object,maxit=100,tol=1e-5,verbose=FALSE,...) {
  if(!is(object,"hmModel")) stop("object must be 'hmModel'")

  # pseudocode
  ns <- object@nstates

  LL <- logLik(object)
  
  converge <- FALSE
  j <- 0

  A <- object@trans
  while(j <= maxit & !converge) {
    for(i in 1:ns) {
        A[,,i] <- predict(object@trModels[[i]])
    }
    B <- exp(apply(object@logdens,c(1,3),sum))
    # TODO: add functionality for inModel
    init <- predict(object@initModel)
    LL.old <- LL
    j <- j+1

    # expectation
    fbo <- fb(init=init,A=A,B=B,ntimes=object@ntimes)
    LL <- fbo$logLike
    
    # maximization
    
    #object@init <- fit(object@init,ip=fbo$gamma[1,])
    #object@init <- matrix(fbo$gamma[1,],nrow=1)
		

    object@initModel <- setpars(object@initModel,values=object@initModel@family$linkfun(fbo$gamma[1,],base=object@initModel@family$base))
    # This should become:
    # lt <- length(object@ntimes)
		# et <- cumsum(object@ntimes)
		# bt <- c(1,et[-lt]+1)
    # object@initModel@y <- fbo$gamma[bt,]
    # object@initModel <- fit(object@initModel)
    
    for(i in 1:ns) {
      object@trModels[[i]]@y <- fbo$xi[,,i]/fbo$gamma[,i]
      object@trModels[[i]] <- fit(object@trModels[[i]],w=as.matrix(fbo$gamma[,i]),ntimes=object@ntimes) # check this
      #object@trModels[[i]] <- fit(object@trModels[[i]],w=NULL,ntimes=object@ntimes) # check this
      #object@trans[,,i] <- exp(logDens(object@trModels[[i]]))
      object@trans[,,i] <- predict(object@trModels[[i]])
      for(k in 1:object@nresp) {
        object@rModels[[i]][[k]] <- fit(object@rModels[[i]][[k]],w=fbo$gamma[,i])
        object@logdens[,k,i] <- logDens(object@rModels[[i]][[k]])
      }
    }
    #object <- setpars(object,getpars(object)) # set parameters and recompute (bit of a roundabout way)
    LL <- logLik(object)
    if(verbose) cat("iteration",j,"logLik:",LL,"\n")
    if( (LL >= LL.old) & (LL - LL.old < tol))  converge <- TRUE
  }
  object
}
