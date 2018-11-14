#
# Ingmar Visser, 8-11-2019
#
# Description
# 
# Computes confidence intervals for (dep)mix model parameters. 
# 
# Details
# 
# Confidence intervals are computed through the variance-covariance matrix
# which in turn is computed using the hessian and the linear constraints
# of the model. See ?vcov, ?hessian and ?standarderror for more details on
# these underlying functions. The confidence intervals are computed using 
# the normal theory assumptions. The desired significance level can be  
# supplied through the 'level' argument. 
# 
# Value
# 
# 

setMethod("confint", "mix",
    function(object, level=0.95, digits=4, fixed=NULL, equal=NULL, 
	conrows=NULL, conrows.upper=NULL, conrows.lower=NULL, 
	tolerance=1e-9, 	
	method="finiteDifferences", ...) {
		
	vc <- vcov(object,fixed=fixed,equal=equal,
		conrows=conrows,conrows.upper=conrows.upper,conrows.lower=conrows.lower,
		tolerance=tolerance,method=method, ...)
		
	ses <- sqrt(diag(vc$vcov))
	
	pars <- getpars(object)
	
	elements <- vc$elements
	
	parsinc <- pars[which(elements=="inc")]
	
	upper <- parsinc+qnorm(0.5+level/2)*ses
	lower <- parsinc-qnorm(0.5+level/2)*ses
	
	ret <- data.frame(pars=round(pars,digits), constr=elements, lower=NA, upper=NA)
		
	ret$lower[which(elements=="inc")] <- round(lower, digits)
	ret$upper[which(elements=="inc")] <- round(upper, digits)
		
	colnames(ret)[3:4] <- c(format.perc(0.5-level/2, 3),format.perc(0.5+level/2, 3))
	
	return(ret)
}
)
