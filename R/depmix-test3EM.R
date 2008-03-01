setwd("/Users/ivisser/Documents/projects/depmixProject/depmixNew/code/depmix/trunk/")

source("depmixS4.R")
source("classes.R")
source("hmModel.R")
source("trGLM.r")
source("EM.R")

load("speed.Rda")

rModels <- list(
  list(
	rModel(formula=rt~1,data=speed,family=gaussian(),pstart=c(5.52,.2)),
	rModel(formula=corr~1,data=speed,family=multinomial())),
  list(
	rModel(formula=rt~1,data=speed,family=gaussian(),pstart=c(6.39,.24)),
	rModel(formula=corr~1,data=speed,family=multinomial(),pstart=c(.098,.902)))
)

trstart=c(0.896,0.104,0.084,0.916)
instart=c(.5,.5)

trstart=c(trstart[1:2],0,0.01,trstart[3:4],0,0)

mod <- depmix(rModels=rModels,data=speed,transition=~Pacc,trstart=trstart,instart=instart)

logLik(mod)

mod@trModels[[1]] <- mod@trModels[[2]] <- trGLM(~Pacc,data=speed,nstate=2)
mod@trModels[[1]]@parameters$coefficients[,2] <- c(-2.153550,.01)
mod@trModels[[2]]@parameters$coefficients[,2] <- c(2.389200,0)

object <- mod
maxit=100
tol=1e-5

fmod <- em(mod,verbose=T)