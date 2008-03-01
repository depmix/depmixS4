
require(VGAM)

setClass("rVGLM",
  contains="rModel"
)

setMethod("fit","rVGLM",
  function(object,w,...) {
    pars <- object@parameters
    fit <- vglm(object@y~object@x,object@family,weight=w,...)
    pars$coefficients <- t(fit$coefficients)  # TODO: setpars gets matrix in wrong order!!! Fix this in setpars.
    object <- setpars(object,unlist(pars))
    object
	}
)

setClass("rMultinom_vglm",
  contains="rVGLM",
  representation(
    base="integer"
  )
)

setMethod("fit","rVGLM",
	function(object,w,...) {
    pars <- object@parameters
    base <- object@base
    y <- cbind(object@y[,-base],object@y[,base])
    x <- object@x
    nan <- unique(apply(y,2,function(x) which(is.na(x))))
    if(length(nan)>0) {
  	  y <- y[-nan,]
      x <- x[-nan,]
	  }
    if(!is.null(w)) {
	    if(length(nan) > 0) w <- w[-nan,]
	    #w <- w[,-base] # TODO: adjust this for more than two levels
      #fit <- vglm(y~x-1,object@family,data=data.frame(y=y,x=x),weight=w,...) 
      fit <- vglm(formula=y~x-1,family=object@family,weight=w) # need to evaluate ...! 
	  } else fit <- vglm(y~x-1,object@family,...)
    pars$coefficients[,-base] <- t(matrix(fit@coefficients,ncol=ncol(x)))  # TODO: setpars gets matrix in wrong order!!! Fix this in setpars.
    object <- setpars(object,unlist(pars))
    object
	}
)


setMethod("predict","rVGLM",
  function(object) {
    pars <- object@parameters$coefficients[,-object@base]
    out <- object@family@inverse(object@x%*%pars)
    tmp <- out[,object@base]
    out[,object@base] <- out[,ncol(out)]
    out[,ncol(out)] <- tmp
    out
  }
)

setMethod("setpars","rVGLM",
  function(object,values,which="all") {
  	npar <- npar(object)
  	if(length(values)!=npar) stop("length of 'values' must be",npar)
  	if(object@family@vfamily[1]=="multinomial") {
  		object@parameters$coefficients[,1] <- values[1:ncol(object@parameters$coefficients)]
		  values <- matrix(values,,ncol(object@x),byrow=TRUE)
  		if(ncol(object@x)>1) object@parameters$coefficients[,2:ncol(object@x)] <- values[2:ncol(object@x),]
  	} else {
  		object@parameters$coefficients <- values[1:length(object@parameters$coefficients)]
  	}
  	if(length(unlist(object@parameters))>length(object@parameters$coefficients)) {
  		if(object@family@vfamily=="gaussian") object@parameters$sd <- values[(length(object@parameters$coefficients)+1)]
  	}
  	return(object)
  }
)

setMethod("getpars","rVGLM",
  function(object,which="par") {
  	parameters <- numeric()
  	if(object@family@vfamily[1]=="multinomial") {
  		# coefficient is usually a matrix here 		
  		parameters <- c(t(object@parameters$coefficients)) # Why transpose?
  	} else {
  		parameters <- unlist(object@parameters)
  	}
  	return(parameters)
  }
)


setClass("trMultinom_vglm",
  contains=c("rMultinom_vglm","trinModel")
)

setMethod("logDens","trMultinom_vglm",
	function(object) {
    log(predict(object))
	}
)

vmultinomial <- function() {
  new("vglmff",
    blurb = "Multinomial logit model\n\n Links:    log(mu[,j]/mu[,M+1])\n Variance: mu[,j]*(1-mu[,j]); -mu[,j]*mu[,k]",
    constraints = expression({
      constraints = cm.vgam(matrix(1, M, 1), x, FALSE, constraints,intercept.apply = FALSE)
      constraints = cm.zero.vgam(constraints, x, NULL, M)
      constraints = cm.nointercept.vgam(constraints, x, NULL, M)
    }),
    deviance = function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
      if(ncol(y) == 1 || ncol(mu) == 1) stop("y and mu must have at least 2 columns")
      double.eps = .Machine$double.eps
      devy = y
      nz = y != 0
      #devy[nz] = y[nz] * log(y[nz])
      #devmu = y * log(mu)
      devmu = mu
      if(any(small <- mu * (1 - mu) < double.eps)) {
        warning("fitted values close to 0 or 1")
        smu = mu[small]
        sy = y[small]
        smu = ifelse(smu < double.eps, double.eps, smu)
        #devmu[small] = sy * log(smu)
        devmu[small] = smu
      }
      devi = (devy - devmu)^2/devy
      if(residuals) {
        M = if(is.matrix(eta)) ncol(eta) else 1
        if(M > 1) return(NULL)
        devi = devi %*% rep(1, ncol(devi))
        return(c(sign(y[, 1] - mu[, 1]) * sqrt(abs(devi) * w)))
      } else sum(w * devi)
    },
    initialize=expression({
      delete.zero.colns = TRUE
      eval(process.categorical.data.vgam)
      M = ncol(y) - 1
      predictors.names = paste("log(mu[,", 1:M, "]/mu[,", M + 1, 
          "])", sep = "")
      y.names = paste("mu", 1:(M + 1), sep = "")
    }),
    inverse = function (eta, extra = NULL) {
      if (any(is.na(eta))) 
          warning("there are NAs in eta in slot inverse")
      phat = cbind(exp(eta), 1)
      ans = phat/as.vector(phat %*% rep(1, ncol(phat)))
      if (any(is.na(ans))) 
          warning("there are NAs here in slot inverse")
      ans
    },
    last = expression({
      misc$link = "mlogit"
      misc$earg = list(mlogit = list())
      dy = dimnames(y)
      if (!is.null(dy[[2]])) 
          dimnames(fit$fitted.values) = dy
    }),
    link = function (mu, extra = NULL) {
      log(mu[, -ncol(mu)]/mu[, ncol(mu)])
    },
    loglikelihood = function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
      if(residuals) stop("loglikelihood residuals not implemented yet") else sum(w * y * log(mu))
    },
    vfamily=c("multinomial","vcategorical"),
    deriv = expression({
      w * (y[, -ncol(y)] - mu[, -ncol(y)])
    }),
    weight = expression({
      tiny = (mu < .Machine$double.eps^0.5) | (mu > 1 - .Machine$double.eps^0.5)
      if(M == 1) { 
          wz = mu[, 1] * (1 - mu[, 1])
      } else {
          index = iam(NA, NA, M, both = TRUE, diag = TRUE)
          wz = -mu[, index$row] * mu[, index$col]
          wz[, 1:M] = wz[, 1:M] + mu[, 1:M]
      }
      atiny = (tiny %*% rep(1, ncol(mu))) > 0
      if(any(atiny)) {
        if(M == 1) {
          wz[atiny] = wz[atiny] * (1 + .Machine$double.eps^0.5) + .Machine$double.eps
        } else {
          wz[atiny, 1:M] = wz[atiny, 1:M] * (1 + .Machine$double.eps^0.5) + .Machine$double.eps
        }
      }
      w * wz
    })
  )
}
  


#VGAM::multinomial
# adapted from gaussian()

trVGLM <- function(formula,family=vmultinomial(),data,base=as.integer(1),nstate,...) {
  call <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data"), names(mf), 0)
	mf <- mf[c(1, m)]
	mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  x <- model.matrix(attr(mf, "terms"),mf)
  y <- matrix(1/nstate,ncol=nstate,nrow=nrow(x))
  #y <- model.response(mf)
  if(!is.matrix(y)) y <- matrix(y,ncol=1)
  parameters <- list()
  parameters$coefficients <- matrix(0,nrow=ncol(x),ncol=nstate)
  # use (weighted) Chi-square deviance
  new("trMultinom_vglm",
    base=base,
    formula=formula,
    family=family,
    parameters=parameters,
    npar=length(unlist(parameters)),
    y=y,
    x=x)
  #  form
}


