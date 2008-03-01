
# the main function to call

# reponse is a list of formulae for each of the response variables,
# defining a glm type equation for each of the variables conditional on the
# state, these formulae may be one-sided, ie only having a response and no
# covariates, family is a list of density functions for the response
# variables

# init is a one-sided formula for the covariates on the initial state
# probabilities or class probabilities which are modelled as a multinomial
# logistic, transition is a one-sided formula for the covariates on the
# transition parameters

# data is an optional dataframe to interpret the variables and
# covariates in; ntimes is a list of individual sequence lengths

# ..start are the starting values for the parameters
# conpat is a vector specifying general linear constraints on the parameters
# 
# should return: an optimized model 
# return: the likelihood of the model
# 

depmix <- function(rModels,transition=~1,init=~1,data=NULL,initdata=NULL,
	trstart=NULL,instart=NULL,prob=TRUE,conpat=NULL,ntimes=NULL,base=1,...) {
    # rModels should be a list of lists, dimension rModels[[nstates]][[nresp]], with 
    # each element of class "rModel"
	
	# maybe move all these checks and the construction of the transition models to hmModel as well	
	# check wheter everything is well-formed and all that
	nstates <- length(rModels)
	nresp <- length(rModels[[1]])
	
	if(!all(lapply(unlist(rModels),is,"rModel"))) stop("'rModels' must be of class 'rModel'")
	if(!all(lapply(rModels, length)==nresp)) stop("number of response variables in rModels differs")	
	
	tst <- FALSE
	if(!is.null(trstart)) {
		tst <- TRUE
		trstart <- matrix(trstart,nstates,byrow=TRUE)
	}
	
	# it may also be possible to do something like this for rModels, especially if 
	# they all have the same form, which would usually be the case anyway
	trModel <- list()
	stationary=FALSE
	if(transition==~1) stationary=TRUE
	for(i in 1:nstates) {
		if(tst) {
			if(stationary) trModel[[i]] <- trinModel(transition,multinomial(base=base),data=data[1,,drop=FALSE],nstates,pstart=trstart[i,],prob=prob)
			else trModel[[i]] <- trinModel(transition,multinomial(base=base),data=data,nstates,pstart=trstart[i,],prob=prob)
		} else {
			if(stationary) trModel[[i]] <- trinModel(transition,multinomial(base=base),data=data[1,,drop=FALSE],nstates,prob=FALSE)
			else trModel[[i]] <- trinModel(transition,multinomial(base=base),data=data,nstates,prob=FALSE)
		}
	}
	
	if(is.null(attributes(data)$ntimes)&is.null(ntimes)) {
		ntimes <- nrow(data)
	} else {
		if(is.null(ntimes)) ntimes <- attributes(data)$ntimes
	}
	
	# initial probabilities model, depending on covariates init(=~1 by default)
	if(init==~1) {
		if(is.null(instart)) {
			initModel <- trinModel(init,data=data.frame(rep(1,length(ntimes))),nst=nstates,family=multinomial())
		} else {
			initModel <- trinModel(init,data=data.frame(rep(1,length(ntimes))),nst=nstates,family=multinomial(),pstart=instart)
		}
	} else {
		if(is.null(initdata)) {
			stop("'Argument initdata missing while the init model is non-trivial")
		} else {
			if(is.null(instart)) {
				initModel <- trinModel(init,data=initdata,nst=nstates,family=multinomial())
			} else {
				initModel <- trinModel(init,data=initdata,nst=nstates,family=multinomial(),pstart=instart)
			}
		}	
	}
	
	mod <- hmModel(rModels,trModel,initModel,ntimes,STATION=stationary)
	
	return(mod)
}
