
# 
# Started by Ingmar Visser, 27-2-2008
# 

# 
# CLASSES & METHODS FOR DEPMIX OBJECTS AND FIT FUNCTIONS
# 

# 
# DEPMIX CLASS
# 

# What does the user want to know after having defined a model?

# 1) the models: formulae or otherwise and the parameter values
# 2) whether parameters are feasible: ie does logLik exist or return nan (very likely with bad parameters)

# What do we need to be able to optimize the model?
# 1) density of the responses
# 2) density of the transitions
# 3) density of the priors
# 4) ntimes: the lengths of individual time series
# 5) whether the models are stationary, ie time-dependent parameters or not

setClass("depmix",
	representation(response="list",
		transition="list",
		# these are the initial state probabilities or priors as they are called in eg mixture or lca models
		prior="ANY",
		dens="array", # B 
		trDens="array", # A
		init="array", # usually called pi 
		stationary="logical",
		ntimes="numeric",
		nstates="numeric",
		nresp="numeric",
		npars="numeric" # number of parameters
	)
)

# 
# METHODS
# 

# Which methods are needed for this class? 
# 1) construction methods: by inputting various models for responses, transitions and priors
# 2) accessor methods: getpars, getdf, nlin, npars, nobs (etc)
# 3) setpars method
# 4) gof methods: logLik, AIC, BIC (etc) 
# 5) print method
# 6) summary method: I am not quite sure what the differences are between print and summary methods?


# CONSTRUCTORS

# the main function constructing a depmix model with full information, ie all models already in place
# this function is probably not ever called by users

setGeneric("depmix", function(response=any, transition=any, ...) standardGeneric("depmix"))

setMethod("depmix",
	signature(response = "list", transition= "list"),
	function(response, transition, prior, ntimes=NULL, stationary=TRUE, ...) {
		
		nstates <- length(response)
		nresp <- length(response[[1]])
		
		# make appropriate ntimes
		if(is.null(ntimes)) {
			ntimes <- nrow(response[[1]][[1]]@y)
		}
		
		# count the number of parameters	
		npars <- npar(prior) 
		for(i in 1:nstates) {
			npars <- npars + sum(sapply(response[[i]],npar))
		}
		npars <- npars + sum(sapply(transition,npar))
				
		# make appropriate array for transition densities
		nt <- sum(ntimes)
		if(stationary) trDens <- array(0,c(1,nstates,nstates))
		else trDens <- array(0,c(nt,nstates,nstates))
		
		# make appropriate array for response densities
		dens <- array(,c(nt,nresp,nstates))
				
		# compute observation and transition densities
		for(i in 1:nstates) {
			for(j in 1:nresp) {
				dens[,j,i] <- dens(response[[i]][[j]]) # remove this response as an argument from the call to setpars
			}
			trDens[,,i] <- dens(transition[[i]])
		}
		
		# compute initial state probabilties
		init <- dens(prior)
		
		new("depmix",response=response,transition=transition,prior=prior,
			dens=dens,trDens=trDens,init=init,stationary=stationary,
			ntimes=ntimes,nstates=nstates,nresp=nresp,npars=npars)
		
	}
)

#
# UNIVARIATE AND MULTIVARIATE MARKOV MIXTURE OF GLM'S
# 

setMethod("depmix",
	signature(response="ANY",transition="formula"),
	function(response,data=NULL,nstates,transition,family=gaussian(),prior=~1,initdata=NULL,
		respstart=NULL,trstart=NULL,instart=NULL,ntimes=NULL, ...) {
		
		if(is.null(data)) {
			if(is.null(ntimes)) stop("'ntimes' must be provided if not in the data")
		} else {
			if(is.null(attr(data,"ntimes"))) {
				if(is.null(ntimes)) ntimes <- nrow(data)
			} else {
				ntimes <- attr(data,"ntimes")
			}
			if(sum(ntimes)!=nrow(data)) stop("'ntimes' and data do not match")
		}
		
		# make response models
		response <- makeResponseModels(response=response,data=data,nstates=nstates,family=family,values=respstart)
		
		# make transition models
		stationary=FALSE
		if(transition==~1) stationary=TRUE
		transition <- makeTransModels(nstates=nstates,formula=transition,data=data,stationary=stationary,values=trstart)
		
		# make prior model
		prior <- makePriorModel(nstates=nstates,ncases=length(ntimes),formula=prior,data=initdata,values=instart)
		
		# call main depmix with all these models, ntimes and stationary
		model <- depmix(response=response,transition=transition,prior=prior,ntimes=ntimes,stationary=stationary)
		
		# deal with starting values here!!!!!!
		
		return(model)
	}
)

setMethod("depmix",
	signature(response="ANY",transition="missing"),
	function(response,data=NULL,nstates,transition=~1,family=gaussian(),prior=~1,initdata=NULL,
		respstart=NULL,trstart=NULL,instart=NULL,ntimes=NULL, ...) {
		
		model <- depmix(response=response,data=data,nstates=nstates,transition=~1,family=family,
			prior=prior,initdata=initdata,respstart=respstart,trstart=trstart,instart=instart,ntimes=ntimes, ...)
		
		return(model)
		
	}
)

makeResponseModels <- function(response,data=NULL,nstates,family,values=NULL,...) {
	
	resp <- response
	response <- list()
	
	st=FALSE
	if(!is.null(values)) st=TRUE
	
	# univariate response data
	if(class(resp)=="formula") {
		nresp <- 1
		for(i in 1:nstates) {
			response[[i]] <- list()
			response[[i]][[1]] <- GLMresponse(resp,data=data,family=family)
			if(st) {
				bp <- npar(response[[i]][[1]])
				response[[i]][[1]] <- GLMresponse(resp,data=data,family=family,pstart=values[1:bp])
				bp <- bp+1
				values <- values[bp:length(values)]
			}
		}
	}
	
	# multi variate response data
	if(is.list(resp)) {
		nresp <- length(resp)
		for(i in 1:nstates) {
			response[[i]] <- list()
			for(j in 1:nresp) {
				response[[i]][[j]] <- GLMresponse(resp[[j]],data=data,family=family[[j]])
				if(st) {
					bp <- npar(response[[i]][[j]])
					response[[i]][[j]] <- GLMresponse(resp[[j]],data=data,family=family[[j]],pstart=values[1:bp])
					bp <- bp+1
					values <- values[bp:length(values)]
				}
			}
		}
	}
	
	return(response)
}

makeTransModels <- function(nstates,formula=~1,data=NULL,stationary,values=NULL, ...) {
	
	# defaults that possibly need some work at some point FIX ME
	base=1
	prob=TRUE
		
	if(!stationary&is.null(data)) stop("non-stationary transition models needs data argument")
		
	# starting values	
	tst <- FALSE
	if(!is.null(values)) {
		tst <- TRUE
		values <- matrix(values,nstates,byrow=TRUE)
	}
		
	models <- list()
	for(i in 1:nstates) {
		if(tst) {
			# in the case of stationary models we do not need data?!?!?!?!?!
			if(stationary) models[[i]] <- transInit(formula,multinomial(base=base),data=data.frame(1),nstates=nstates,pstart=values[i,],prob=prob)
			else models[[i]] <- transInit(formula,multinomial(base=base),data=data,nstates=nstates,pstart=values[i,],prob=prob)
		} else {
			if(stationary) models[[i]] <- transInit(formula,multinomial(base=base),data=data.frame(1),nstates=nstates,prob=FALSE)
			else models[[i]] <- transInit(formula,multinomial(base=base),data=data,nstates=nstates,prob=FALSE)
		}
	}
	
	return(models)
}

makePriorModel <- function(nstates,ncases,formula=~1,data=NULL,values=NULL) {
	
	# these arguments need to be added at some point FIX ME
	base=1
	prob=TRUE
	
	# initial probabilities model, depending on covariates init(=~1 by default)
	if(formula==~1) {
		initModel <- transInit(~1,data=data.frame(rep(1,ncases)),nst=nstates,family=multinomial(),pstart=values)
	} else {
		if(is.null(data)) {
			stop("'Argument initdata missing while the init model is non-trivial")
		} else {
			initModel <- transInit(formula,data=data,nst=nstates,family=multinomial(),pstart=values)
		}	
	}
	
	return(initModel)
}

setMethod("show","depmix",
	function(object) {
		cat("Initial state probabilties model \n")
		print(object@prior)
		cat("\n")
		for(i in 1:object@nstates) {
			cat("Transition model for state (component)", i,"\n")
			print(object@transition[[i]])
			cat("\n")
		}
		cat("\n")
		for(i in 1:object@nstates) {
			cat("Response model(s) for state", i,"\n\n")
			for(j in 1:object@nresp) {
				cat("Response model for response",j,"\n")
				print(object@response[[i]][[j]])
				cat("\n")
			}
			cat("\n")
		}
	}
)

setMethod("getpars","depmix",
	function(object,which="pars",...) {
		parameters <- getpars(object@prior,which=which)
		for(i in 1:object@nstates) {
			parameters <- c(parameters,getpars(object@transition[[i]],which=which))
		}
		for(i in 1:object@nstates) {
			for(j in 1:object@nresp) {
				parameters <- c(parameters,getpars(object@response[[i]][[j]],which=which))
			}
		}
		return(parameters)
	}
)

setGeneric("nobs", function(object, ...) standardGeneric("nobs"))

setMethod("nobs", signature(object="depmix"),
	function(object, ...) {
		sum(object@ntimes)
	}
)

# accessor functions
setMethod("npar","depmix",
	function(object) return(object@npars)
)

setGeneric("freepars", function(object, ...) standardGeneric("freepars"))

# depends on nlin(object) and getpars(object)
setMethod("freepars","depmix",
	function(object) {
		free <- sum(!getpars(object,which="fixed"))
# 		free <- free-nlin(object)
		free
	}
)

setGeneric("logLik", function(object, ...) standardGeneric("logLik"))

source("lystig.R")
source("fb.R")

# depends on getpars and nobs
setMethod("logLik",signature(object="depmix"),
	function(object,method="lystig") { # TODO: initial dens and response densities are recomputed here, but this is also done in setpars at least for the response densities !!!!!!!!
		if(method=="fb") ll <- fb(object@init,object@trDens,apply(object@dens,c(1,3),prod),object@ntimes,object@stationary)$logLike
		if(method=="lystig") ll <- lystig(object@init,object@trDens,apply(object@dens,c(1,3),prod),object@ntimes,object@stationary)$logLike
		attr(ll, "df") <- freepars(object)
		attr(ll, "nobs") <- nobs(object)
		class(ll) <- "logLik"
		ll
	}
)

# depends on logLik and freepars
setMethod("AIC", signature(object="depmix"),
	function(object, ..., k=2){
		c(-2 * logLik(object) + freepars(object) * k)
	}
)

setGeneric("BIC", function(object, ...) standardGeneric("BIC"))

# depends on logLik, freepars and nobs
setMethod("BIC", signature(object="depmix"),
	function(object, ...){
		c(-2 * logLik(object) + freepars(object) * log(nobs(object)))
	}
)

# depends on npar
setMethod("setpars","depmix",
	function(object,values,which="pars",...) {
		if(!(length(values)==npar(object))) stop("Argument 'values' has incorrect length")
		bp <- npar(object@prior)
		switch(which,
			"pars" = {
				if(!all(getpars(object@prior,which=which) == values[1:bp])) {
					object@prior=setpars(object@prior,values[1:bp],which=which)
					# recompute init probabilities
					object@init <- dens(object@prior)
				}
			},
			"fixed" = {
				object@prior <- setpars(object@prior,values[1:bp],which=which)
			}
		)
		bp <- bp+1
		values <- values[bp:npar(object)]
		for(i in 1:object@nstates) {
			bp <- npar(object@transition[[i]])
			switch(which,
				"pars"= {
					if(!all(getpars(object@transition[[i]]) == values[1:bp])) {
						object@transition[[i]] <- setpars(object@transition[[i]],values[1:bp])
						# recompute transition densities if pars have changed
						object@trDens[,,i] <- dens(object@transition[[i]])
					}
				},
				"fixed" = {
					object@transition[[i]] <- setpars(object@transition[[i]],values[1:bp],which="fixed")
				}
			)
			bp <- bp+1
			values <- values[bp:length(values)]
		}
		for(i in 1:object@nstates) {
			for(j in 1:object@nresp) {
				bp <- npar(object@response[[i]][[j]])
				switch(which,
					"pars" = {
						if(!all(getpars(object@response[[i]][[j]]) == values[1:bp])) {
							object@response[[i]][[j]] <- setpars(object@response[[i]][[j]],values[1:bp])
							# recompute observation densities if pars have changed
							object@dens[,j,i] <- dens(object@response[[i]][[j]])
						}
					},
					"fixed" = {
						object@response[[i]][[j]] <- setpars(object@response[[i]][[j]],values[1:bp],which="fixed")
					}
				)	
				bp <- bp+1
				values <- values[bp:length(values)]
			}
		}			
		return(object)
	}
)



