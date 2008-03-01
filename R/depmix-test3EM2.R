setwd("/Users/ivisser/Documents/projects/depmixProject/depmixNew/code/depmix/trunk/")

source("depmixS4.R")
source("classes.R")
source("hmModel.R")
source("lystig.R")
source("fb.R")
source("trGLM.r")
source("rVGLM.r")

source("EM.R")

maxit=100
tol=1e-5

load("speed.Rda")

rModels <- list(
  list(
	rModel(formula=rt~1,data=speed,family=gaussian(),pstart=c(5.5,.2))),
  list(
	rModel(formula=rt~1,data=speed,family=gaussian(),pstart=c(6.3,.2)))
)

trstart=c(0.8,0.2,0.1,0.9)
instart=c(.5,.5)


mod <- depmix(rModels=rModels,data=speed,transition=~1,trstart=trstart,instart=instart)

logLik(mod)

fmod <- em(mod,verbose=T)

# 
# OPTIMIZE USING RDONLP
# 

fixed <- getpars(mod,"fixed")
allpars <- getpars(mod)

pars <- allpars[!fixed]

logl <- function(pars) {
	allpars[!fixed] <- pars
	mod <- setpars(mod,allpars)
	-logLik(mod)
}

library(Rdonlp2)

cntrl <- donlp2.control(hessian=FALSE,difftype=1,epsfcn=1e-6)

res1 <- donlp2(pars,logl,control=cntrl)

allpars[!fixed] <- res1$par
mod <- setpars(mod,allpars)

ll <- logLik(mod)


# final loglike
# [1] -84.3424

# pars:   
# init: 	0.00        20.58             
# trans: 	0.00 		-2.13            
# 			0.00        2.39
#         
# obser:	5.51     	0.19    
# 			6.38        0.24




		