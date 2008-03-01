

setClass("trGLM",
  contains="rGLM",
  representation(
    base="numeric"
  )
)

trGLM <- function(formula,family=binomial(),data,base=as.integer(1),nstate,...) {
  call <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data"), names(mf), 0)
	mf <- mf[c(1, m)]
	mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  x <- model.matrix(attr(mf, "terms"),mf)
  y <- matrix(1/nstate,ncol=nstate,nrow=nrow(x))
  if(!is.matrix(y)) y <- matrix(y,ncol=1)
  parameters <- list()
  parameters$coefficients <- matrix(0,nrow=ncol(x),ncol=nstate)
  new("trGLM",
    base=base,
    formula=formula,
    family=family,
    parameters=parameters,
    npar=length(unlist(parameters)),
    y=y,
    x=x)
  #  form
}

setMethod("setpars","trGLM",
  function(object,values,which="all") {
  	npar <- npar(object)
  	if(length(values)!=npar) stop("length of 'values' must be",npar)
#  	if(object@family@vfamily[1]=="multinomial") {
  		object@parameters$coefficients[,1] <- values[1:ncol(object@parameters$coefficients)]
		  values <- matrix(values,,ncol(object@x),byrow=TRUE)
  		if(ncol(object@x)>1) object@parameters$coefficients[,2:ncol(object@x)] <- values[2:ncol(object@x),]
#  	} else {
#  		object@parameters$coefficients <- values[1:length(object@parameters$coefficients)]
#  	}
  	if(length(unlist(object@parameters))>length(object@parameters$coefficients)) {
  		if(object@family@vfamily=="gaussian") object@parameters$sd <- values[(length(object@parameters$coefficients)+1)]
  	}
  	return(object)
  }
)

setMethod("getpars","trGLM",
  function(object,which="par") {
  	parameters <- numeric()
#  	if(object@family@vfamily[1]=="multinomial") {
  		# coefficient is usually a matrix here
  		parameters <- c(t(object@parameters$coefficients)) # Why transpose?
#  	} else {
 # 		parameters <- unlist(object@parameters)
 # 	}
  	return(parameters)
  }
)

setMethod("fit","trGLM",
	function(object,w,...) {
    pars <- object@parameters
    base <- object@base
    y <- as.matrix(object@y[,-base])
    x <- as.matrix(object@x)
    nan <- unique(apply(y,2,function(x) which(is.na(x))))
    if(length(nan)>0) {
  	  y <- y[-nan,]
      x <- x[-nan,]
	  }
    if(!is.null(w)) {
	    if(length(nan) > 0) w <- w[-nan,]
	    #w <- w[,-base] # TODO: adjust this for more than two levels
      #fit <- vglm(y~x-1,object@family,data=data.frame(y=y,x=x),weight=w,...)
      fit <- glm(formula=y~x-1,family=object@family,weight=w) # need to evaluate ...!
	  } else fit <- glm(y~x-1,object@family)
    pars$coefficients[,-base] <- fit$coefficients  # TODO: setpars gets matrix in wrong order!!! Fix this in setpars.
    #object <- setpars(object,unlist(pars))
    object@parameters <- pars
    object
	}
)

setMethod("predict","trGLM",
  function(object) {
    out <- matrix(0,nrow=nrow(object@y),ncol=ncol(object@y))
    pars <- object@parameters$coefficients[,-object@base]
    out[,-object@base] <- object@family$linkinv(object@x%*%pars)
    out[,object@base] <- 1-rowSums(out)
#    tmp <- out[object@base]
#    out[,object@base] <- out[,ncol(out)]
#    out[,ncol(out)] <- tmp
    out
  }
)