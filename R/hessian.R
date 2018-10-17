# Ingmar Visser, 17 oktober 2018
# 
# Description
# The function hessian computes the hessian of a (dep)mix object at the current 
# parameter values; it has a method argument; the only method currently implemented
# is the finite differences method. 
#
# Details
# 
# The function optionally accepts arguments related to linear constraints, which are 
# neccessary to have when the hessian is used to compute the variance-covariance 
# matrix of the parameters. 
# 
# The function also checks whether parameters are estimated on the boundary and 
# leaves them out of the process when this is the case. 
# 
# Fixed parameters are similarly ignored when computing the hessian
# 
# Value
#
# The function returns a named list with the following elements:
# 
# hessian: the hessian of the parameters
# 
# elements: vector of length npar(object) indicating for each parameter whether it 
# is 'inc'luded, 'fix'ed, or estimated on the boundary, 'bnd'; the dimension of the hessian 
# is thus the number of non-fixed parameters minus the number of boundary parameters. 
# 
# lincon: the linear constraint matrix needed to compute the variance-covariance 
# matrix; it only contains the parts of the linear constraint matrix that relate to equality 
# constraints; moreover, the columns related to 'fix'ed and boundary ('bnd') parameters
# are left out. 
#

setMethod("hessian", "mix",
    function(object, fixed=NULL, equal=NULL, 
	conrows=NULL, conrows.upper=NULL, conrows.lower=NULL, 
	tolerance=1e-9, 	
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
	
	# incorporate general linear constraints, if any, in lincon matrix
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
	
	# return vector with specification of 'inc'luded, 'fix'ed and 'bnd'ary parameters
	elements <- rep("inc",npar(object))
	
	# identify parameters that are on their boundary
	low <- which(sapply(as.numeric(allpars-par.l),all.equal,tolerance=tolerance,0)==TRUE)
	up <- which(sapply(as.numeric(allpars-par.u),all.equal,tolerance=tolerance,0)==TRUE)
	bnd <- union(low, up)
	
	
	if(length(which(fixed)>0)) elements[which(fixed)] <- "fix"
	if(length(bnd)>0) elements[bnd] <- "bnd"
	
	# get the reduced set of parameters, ie the ones that the hessian will be computed for
	# only non-fixed parameters
	pars <- allpars[which(elements=="inc")]
	
	# select only those columns of the constraint matrix that correspond to non-fixed parameters
	lincon <- lincon[,which(elements=="inc"),drop=FALSE]
	
	# remove redundant rows in lincon (all zeroes)
	allzero <- which(apply(lincon,1,function(y) all(y==0)))
	if(length(allzero)>0) {
		lincon <- lincon[-allzero,,drop=FALSE]
		lin.u <- lin.u[-allzero]
		lin.l <- lin.l[-allzero]
	}
	
	# remove rows of lincon with inequality constraints
	dflu <- lin.u-lin.l
	ineq <- which(dflu!=0)
	if(length(ineq)>0) {
		lincon <- lincon[-ineq,,drop=FALSE]
		lin.u <- lin.u[-ineq]
		lin.l <- lin.l[-ineq]
	}
	
	# make loglike function that only depends on pars
	logl <- function(pars) {
		allpars[which(elements=="inc")] <- pars
		object <- setpars(object,allpars)
		ans = -as.numeric(logLik(object))
		if(is.na(ans)) ans = 1000000 # remove magic number here!!!
		ans
	}
	
	fdh <- fdHess(pars,logl)
	
	# also return list of length npar that specifies for which parameters 
	# the hessian has been computed and for which this has been skipped due
	# to being 1) fixed or 2) on the boundary
		
	return(list(hessian=fdh$Hessian,elements=elements,lincon=lincon))
}
)