
# 
# Ingmar Visser, 11-6-2008
# 

# Changes
# - added lin.upper and lin.lower slots to these objects

# 
# MIX.FITTED CLASS
# 

setClass("mix.fitted",
	representation(message="character", # convergence information
		conMat="matrix", # constraint matrix on the parameters for general linear constraints
		lin.upper="numeric", # upper bounds for linear constraint
		lin.lower="numeric", # lower bounds for linear constraints
		posterior="data.frame" # posterior probabilities for the states
	),
	contains="mix"
)

# accessor functions

setMethod("posterior","mix.fitted",
	function(object) {
		return(object@posterior)
	}
)

setMethod("show","mix.fitted",
	function(object) {
		cat("Convergence info:",object@message,"\n")
		print(logLik(object))
		cat("AIC: ", AIC(object),"\n")
		cat("BIC: ", BIC(object),"\n")
	}
)

setMethod("summary","mix.fitted",
	function(object,which="all") {
		ans=switch(which,
			"all" = 1,
			"response" = 2,
			"prior" = 3,
			stop("Invalid 'which' argument in summary of fitted mix model")
		)
		if(ans==1|ans==3) {
			cat("Mixture probabilities model \n")
			show(object@prior)
			cat("\n")
		}
		if(ans==1|ans==2) {
			for(i in 1:object@nstates) {
				cat("Response model(s) for state", i,"\n\n")
				for(j in 1:object@nresp) {
					cat("Response model for response",j,"\n")
					show(object@response[[i]][[j]])
					cat("\n")
				}
				cat("\n")
			}
		}
	}	
)

# 
# Ingmar Visser, 23-3-2008
# 

# 
# DEPMIX.FITTED CLASS
# 

setClass("depmix.fitted",
	representation(message="character", # convergence information
		conMat="matrix", # constraint matrix on the parameters for general linear constraints
		lin.upper="numeric", # upper bounds for linear constraints
		lin.lower="numeric", # lower bounds for linear constraints
		posterior="data.frame" # posterior probabilities for the states
	),
	contains="depmix"
)

# accessor functions

setMethod("posterior","depmix.fitted",
	function(object) {
		return(object@posterior)
	}
)

setMethod("show","depmix.fitted",
	function(object) {
		cat("Convergence info:",object@message,"\n")
		print(logLik(object))
		cat("AIC: ", AIC(object),"\n")
		cat("BIC: ", BIC(object),"\n")
	}
)

# copied from hmmr (and removed there)

setMethod("summary","depmix.fitted",
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
				if(compact & object@prior@formula==~1) {
						pr <- object@prior@parameters$coefficients
						rownames(pr) <- ""
						colnames(pr) <- paste("St",1:ns,sep="")
						cat(pr,"\n")
				} else show(object@prior)
		}
		if(ans==1|ans==4) {
				# show the transition models
				if(compact & object@transition[[1]]@formula==~1) {
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



