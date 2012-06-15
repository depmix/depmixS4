# depends on getpars and nobs
setMethod("logLik",signature(object="depmix"),
	#function(object,method="lystig") { 
	function(object,method="fb") { #4/5/2012: set to fb as this is now in C
		if(method=="fb") ll <- fb(object@init,object@trDens,object@dens,object@ntimes,object@stationary)$logLike
		if(method=="lystig") ll <- lystig(object@init,object@trDens,object@dens,object@ntimes,object@stationary)$logLike
		attr(ll, "df") <- freepars(object)
		attr(ll, "nobs") <- nobs(object)
		class(ll) <- "logLik"
		ll
	}
)

# depends on getpars and nobs
setMethod("logLik",signature(object="mix"),
	#function(object,method="lystig") { 
	function(object,method="fb") { 
		if(method=="fb") ll <- fb(object@init,matrix(0,1,1),object@dens,object@ntimes,TRUE)$logLike
		if(method=="lystig") ll <- lystig(object@init,matrix(0,1,1),object@dens,object@ntimes,TRUE)$logLike
		attr(ll, "df") <- freepars(object)
		attr(ll, "nobs") <- nobs(object)
		class(ll) <- "logLik"
		ll
	}
)
