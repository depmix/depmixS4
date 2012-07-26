
# 
# Ingmar Visser, 11-6-2008
# 

# 
# DEPMIX CLASS BELOW THE MIX CLASS
# 

# 
# Class definition, accessor functions, print and summary methods
# 

# 
# MIX CLASS
# 

setClass("mix",
	representation(response="list", # response models
		prior="ANY", # the prior model (multinomial)
		dens="array", # response densities (B)
		init="array", # usually called pi 
		nstates="numeric",
		nresp="numeric",
		ntimes="numeric",
		npars="numeric" # number of parameters
	)
)

# accessor functions
setMethod("npar","mix",
	function(object) return(object@npars)
)

setMethod("ntimes","mix",
	function(object) return(object@ntimes)
)

setMethod("nstates","mix",
	function(object) return(object@nstates)
)

setMethod("nresp","mix",
	function(object) return(object@nresp)
)

setMethod("is.stationary",signature(object="mix"),
  function(object) {
		return(TRUE)
	}
)

setMethod("getmodel","mix",
		function(object, which="response", state=1, number) {
				ans=switch(
						which,
						"response" = 1,
						"prior" = 2,
						stop("Invalid 'which' argument in getmodel on 'mix' object.")
				)
				if(ans==1) {
						if(is.null(number)) number <- 1
						return(object@response[[state]][[number]])
				}
				if(ans==2) return(object@prior)
		}
)

setMethod("simulate",signature(object="mix"),
	function(object,nsim=1,seed=NULL,...) {
		
		if(!is.null(seed)) set.seed(seed)
		
		ntim <- ntimes(object)
		nt <- sum(ntim)
		bt <- 1:nt
		
		nr <- nresp(object)
		ns <- nstates(object)
		
		# simulate state sequences first, then observations
		
		# random generation is slow when done separately for each t, so first draw
		# variates for all t, and then determine state sequences iteratively
		states <- array(,dim=c(nt,nsim))
		states[bt,] <- simulate(object@prior,n=nsim,is.prior=T)
		sims <- array(,dim=c(nt,ns,nsim))
				
		states <- as.vector(states)
		responses <- list(length=nr)
		#responses <- array(,dim=c(nt,nr,nsim))
		for(i in 1:nr) {
			tmp <- matrix(,nrow=nt*nsim,ncol=NCOL(object@response[[1]][[i]]@y))
			for(j in 1:ns) {
				tmp[states==j,] <- simulate(object@response[[j]][[i]],nsim=nsim)[states==j,]
			}
			responses[[i]] <- tmp
		}
		
		# generate new mix.sim object
		class(object) <- c("mix.sim")
		object@states <- as.matrix(states)
		
		object@prior@x <- as.matrix(apply(object@prior@x,2,rep,nsim))
		for(j in 1:ns) {
			for(i in 1:nr) {
				object@response[[j]][[i]]@y <- as.matrix(responses[[i]])
				object@response[[j]][[i]]@x <- as.matrix(apply(object@response[[j]][[i]]@x,2,rep,nsim))
			}
		}
		object@ntimes <- rep(object@ntimes,nsim)
		
		# make appropriate array for transition densities
		nt <- sum(object@ntimes)
		
		# make appropriate array for response densities
		dns <- array(,c(nt,nr,ns))
		
		# compute observation and transition densities
		for(i in 1:ns) {
			for(j in 1:nr) {
				dns[,j,i] <- dens(object@response[[i]][[j]]) # remove this response as an argument from the call to setpars
			}
		}
		
		# compute initial state probabilties
		object@init <- dens(object@prior)
		object@dens <- dns
		
		return(object)
	}
)

# setMethod("getModel",signature(object="mix"),
# 	function(object,which="response",...) {
# 		res <- switch(which,
# 			"prior"=object@prior,
# 			"response"=object@response)
# 		res
# 	}
# )

# 
# PRINT method
# 

setMethod("show","mix",
		function(object) {
				summary(object)
		}
)

# 
# SUMMARY method
# 

setMethod("summary","mix",
		function(object,which="all",compact=TRUE) {
				ans=switch(which,
						"all" = 1,
						"response" = 2,
						"prior" = 3,
						stop("Invalid 'which' argument in summary of fitted mix model")
				)
				if(ans==1|ans==3) {
						# show the prior models
						cat("Mixture probabilities model \n")
						if(object@prior@formula==~1) {
								pr <- object@prior@parameters$coefficients
								rownames(pr) <- ""
								colnames(pr) <- paste("St",1:ns,sep="")
								print(pr)
						} else show(object@prior)
						cat("\n")
				}
				if(ans==1|ans==2) {
						# show the response models
						if(!compact) {
								for(i in 1:ns) {
										cat("Response model(s) for state", i,"\n\n")
										for(j in 1:object@nresp) {
												cat("Response model for response",j,"\n")
												show(object@response[[i]][[j]])
												cat("\n")
										}
										cat("\n")
								}
						} else {
								cat("Response models \n")
								for(j in 1:object@nresp) {
										cat("Resp",j, ":", object@response[[1]][[j]]@family$family, "\n")
								}
								cat("Response parameters \n")
								pars <- list()
								np <- numeric(object@nresp)
								for(j in 1:object@nresp) {
										np[j] <- npar(object@response[[1]][[j]])
										pars[[j]] <- matrix(,nr=ns,nc=np[j])
								}
								allpars <- matrix(,nr=ns,nc=0)
								nms <- c()
								for(j in 1:object@nresp) {
										for(i in 1:ns) {
												pars[[j]][i,]=getpars(object@response[[i]][[j]])
										}
										nms <- c(nms,paste("Resp",j,1:np[j],sep="."))
										allpars <- cbind(allpars,pars[[j]])					
								}
								rownames(allpars) <- paste("St",1:ns,sep="")
								colnames(allpars) <- nms
								print(allpars)
						}
				}
		}	
)


# 
# Ingmar Visser, 23-3-2008
# 

# 
# Class definition, accessor functions, print and summary methods
# 

# 
# DEPMIX CLASS
# 

setClass("depmix",
	representation(transition="list", # transition models (multinomial logistic)
		trDens="array", # transition densities (A)
		stationary="logical"
	),
	contains="mix"
)


# 
# PRINT method
# 

setMethod("show","depmix",
		function(object) {
				summary(object)
		}
)

setMethod("getmodel","depmix",
		function(object, which="response", state=1, number) {
				ans=switch(
						which,
						"response" = 1,
						"prior" = 2,
						"transition"=3,
						stop("Invalid 'which' argument in getmodel on 'mix' object.")
				)
				if(ans==1) {
						if(is.null(number)) number <- 1
						return(object@response[[state]][[number]])
				}
				if(ans==2) return(object@prior)
				if(ans==3) return(object@transition[[state]])
		}
)

setMethod("is.stationary",signature(object="depmix"),
  function(object) {
		return(object@stationary)
	}
)

setMethod("simulate",signature(object="depmix"),
	function(object,nsim=1,seed=NULL,...) {
		
		if(!is.null(seed)) set.seed(seed)
		
		ntim <- ntimes(object)
		nt <- sum(ntim)
		lt <- length(ntim)
		et <- cumsum(ntim)
		bt <- c(1,et[-lt]+1)
		
		nr <- nresp(object)
		ns <- nstates(object)
		
		# simulate state sequences first, then observations
		
		# random generation is slow when done separately for each t, so first draw
		#   variates for all t, and then determine state sequences iteratively
		states <- array(,dim=c(nt,nsim))
		states[bt,] <- simulate(object@prior,n=nsim,is.prior=T)
		sims <- array(,dim=c(nt,ns,nsim))
		for(i in 1:ns) {
			if(is.stationary(object)) {
				# TODO: this is a temporary fix!!! 
				sims[,i,] <- simulate(object@transition[[i]],nsim=nsim,times=rep(1,nt))
			} else {
				sims[,i,] <- simulate(object@transition[[i]],nsim=nsim)
			}
		}
		# track states
		for(case in 1:lt) {
			for(i in (bt[case]+1):et[case]) {
				states[i,] <- sims[cbind(i,states[i-1,],1:nsim)]
			}
		}
		
		states <- as.vector(states)
		responses <- list(length=nr)
		#responses <- array(,dim=c(nt,nr,nsim))
		for(i in 1:nr) {
			tmp <- matrix(,nrow=nt*nsim,ncol=NCOL(object@response[[1]][[i]]@y))
			for(j in 1:ns) {
				tmp[states==j,] <- simulate(object@response[[j]][[i]],nsim=nsim)[states==j,]
			}
			responses[[i]] <- tmp
		}
		
		# generate new depmix.sim object
		class(object) <- c("depmix.sim")
		object@states <- as.matrix(states)
		
		object@prior@x <- as.matrix(apply(object@prior@x,2,rep,nsim))
		for(j in 1:ns) {
			if(!is.stationary(object)) object@transition[[j]]@x <- as.matrix(apply(object@transition[[j]]@x,2,rep,nsim))
			for(i in 1:nr) {
				object@response[[j]][[i]]@y <- as.matrix(responses[[i]])
				object@response[[j]][[i]]@x <- as.matrix(apply(object@response[[j]][[i]]@x,2,rep,nsim))
			}
		}
		object@ntimes <- rep(object@ntimes,nsim)
		
		# make appropriate array for transition densities
		nt <- sum(object@ntimes)
		if(is.stationary(object)) trDens <- array(0,c(1,ns,ns)) else trDens <- array(0,c(nt,ns,ns))
		
		# make appropriate array for response densities
		dns <- array(,c(nt,nr,ns))
		
		# compute observation and transition densities
		for(i in 1:ns) {
			for(j in 1:nr) {
				dns[,j,i] <- dens(object@response[[i]][[j]]) # remove this response as an argument from the call to setpars
			}
			trDens[,,i] <- dens(object@transition[[i]])
		}
		
		# compute initial state probabilties
		object@init <- dens(object@prior)
		object@trDens <- trDens
		object@dens <- dns
		
		return(object)
	}
)

# setMethod("getModel",signature(object="depmix"),
# 	function(object,which="response",...) {
# 		res <- switch(which,
# 			"prior"=object@prior,
# 			"response"=object@response,
# 			"transition"=object@transition)
# 		res
# 	}
# )

# 
# SUMMARY method: to do
# 

# copied from hmmr (and removed there)

setMethod("summary","depmix",
	function(object,which="all", compact=TRUE) {
		ns <- object@nstates
		ans=switch(which,
			"all" = 1,
			"response" = 2,
			"prior" = 3,
			"transition" = 4,
			stop("Invalid 'which' argument in summary of fitted depmix model")
		)
		if(ans==1|ans==3) {
				# show the prior models
				cat("Initial state probabilties model \n")
				if(object@prior@formula==~1) {
						pr <- object@prior@parameters$coefficients
						rownames(pr) <- ""
						colnames(pr) <- paste("St",1:ns,sep="")
						print(pr)
						cat("\n")
				} else show(object@prior)
		}
		if(ans==1|ans==4) {
				# show the transition models
				if(object@transition[[1]]@formula==~1) {
						cat("\nTransition matrix \n")
						pars <- getpars(object)
						trm <- matrix(pars[(ns+1):(ns^2+ns)],ns,ns,byr=T)
						rownames(trm) <- paste("fromS",1:ns,sep="")
						colnames(trm) <- paste("toS",1:ns,sep="")
						print(trm)
						cat("\n")
				} else {
						for(i in 1:ns) {
								cat("Transition model for state (component)", i,"\n")
								show(object@transition[[i]])
								cat("\n")
						}
						cat("\n")
				}
		}
		if(ans==1|ans==2) {
				# show the response models
			if(!compact) {
				for(i in 1:ns) {
					cat("Response model(s) for state", i,"\n\n")
					for(j in 1:object@nresp) {
						cat("Response model for response",j,"\n")
						show(object@response[[i]][[j]])
						cat("\n")
					}
					cat("\n")
				}
			} else {
				cat("Response parameters \n")
				for(j in 1:object@nresp) {
						cat("Resp",j, ":", object@response[[1]][[j]]@family$family, "\n")
				}
				pars <- list()
				np <- numeric(object@nresp)
				for(j in 1:object@nresp) {
					np[j] <- npar(object@response[[1]][[j]])
					pars[[j]] <- matrix(,nr=ns,nc=np[j])
				}
				allpars <- matrix(,nr=ns,nc=0)
				nms <- c()
				for(j in 1:object@nresp) {
					for(i in 1:ns) {
						pars[[j]][i,]=getpars(object@response[[i]][[j]])
				}
				nms <- c(nms,paste("Resp",j,1:np[j],sep="."))
					allpars <- cbind(allpars,pars[[j]])					
				}
				rownames(allpars) <- paste("St",1:ns,sep="")
				colnames(allpars) <- nms
				print(allpars)
			}
		}
	}
)


