
setMethod("hessian", "mix",
    function(object,fixed=NULL,equal=NULL,
	conrows=NULL,conrows.upper=NULL,conrows.lower=NULL,	
	method="finiteDifferences", ...) {

	if(is.nan(logLik(object))) stop("Log likelihood is 'NaN'; cannot compute hessian. ")
	
	# check for presence of constraints
	fi <- !is.null(fixed)
	cr <- !is.null(conrows)
	eq <- !is.null(equal)
	
	constr <- any(c(fi,cr,eq))
	
	# determine which parameters are fixed
	if(fi) {
		if(length(fixed)!=npar(object)) stop("'fixed' does not have correct length")
	} else {
		if(eq) {
			if(length(equal)!=npar(object)) stop("'equal' does not have correct length")
			fixed <- !pa2conr(equal)$free
		} else {
			fixed <- getpars(object,"fixed")
		}
	}
	
	# set those fixed parameters in the appropriate submodels
	object <- setpars(object,fixed,which="fixed")
	
	constraints <- getConstraints(object)
	
	lincon=constraints$lincon
	lin.u=constraints$lin.u
	lin.l=constraints$lin.l
	par.u=constraints$par.u
	par.l=constraints$par.l
	
	# incorporate equality constraints provided with the hessian function, if any
	if(eq) {
		if(length(equal)!=npar(object)) stop("'equal' does not have correct length")
		equal <- pa2conr(equal)$conr
		lincon <- rbind(lincon,equal)
		lin.u <- c(lin.u,rep(0,nrow(equal)))
		lin.l <- c(lin.l,rep(0,nrow(equal)))				
	}
	
	# incorporate general linear constraints, if any
	if(cr) {
		if(ncol(conrows)!=npar(object)) stop("'conrows' does not have the right dimensions")
		lincon <- rbind(lincon,conrows)
		if(any(conrows.upper==0)) {
			lin.u <- c(lin.u,rep(0,nrow(conrows)))
		} else {
			if(length(conrows.upper)!=nrow(conrows)) stop("'conrows.upper does not have correct length")
			lin.u <- c(lin.u,conrows.upper)
		}
		if(any(conrows.lower==0)) {
			lin.l <- c(lin.l,rep(0,nrow(conrows)))
		} else {
			if(length(conrows.lower)!=nrow(conrows)) stop("'conrows.lower does not have correct length")
			lin.l <- c(lin.l,conrows.lower)
		}
	}
	
	# get the full set of parameters
	allpars <- getpars(object)
	
	# get the reduced set of parameters, ie the ones that the hessian will be computed for
	# only non-fixed parameters
	pars <- allpars[!fixed]
	
	# TODO: now also remove parameters that are on their boundary
	
	
	# select only those columns of the constraint matrix that correspond to non-fixed parameters
	linconFull <- lincon
	lincon <- lincon[,!fixed,drop=FALSE]
	
	# remove redundant rows in lincon (all zeroes)
	allzero <- which(apply(lincon,1,function(y) all(y==0)))
	if(length(allzero)>0) {
		lincon <- lincon[-allzero,,drop=FALSE]
		lin.u <- lin.u[-allzero]
		lin.l <- lin.l[-allzero]
	}
	
	# TODO: remove rows of lincon with inequality constraints!!!!
	
	
	# make loglike function that only depends on pars
	logl <- function(pars) {
		allpars[!fixed] <- pars
		object <- setpars(object,allpars)
		ans = -as.numeric(logLik(object))
		if(is.na(ans)) ans = 100000 # remove magic number here!!!!!!!!
		ans
	}
	
	fdh <- fdHess(pars,logl)
	
	# also return list of length npar that specifies for which parameters 
	# the hessian has been computed and for which this has been skipped due
	# to being 1) fixed or 2) on the boundary
		
	return(list(hessian=fdh$Hessian,lincon=lincon))
}
)