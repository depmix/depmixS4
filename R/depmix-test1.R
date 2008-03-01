
# 
# Started by Ingmar Visser & Maarten Speekenbrink, 31-1-2007
# 
# Usage: go to trunk directory and source this file in R, if the program
# still works it should return TRUE at every test (or make immediate sense
# otherwise)

# Changes: 
# 
# 16-5-2007, made this a new copy with basic tests of computing the
# likelihood, copying parameters etc
# 
# Other tests with optimization of models are moved to depmix-test2.R
# 

# 
# get original depmix values
# 

setwd("/Users/ivisser/Documents/projects/depmixProject/depmixNew/code/depmix/trunk/")

.libPaths(new="/Users/ivisser/Library/R/library/")
library(depmix)

data(speed)
mod <- dmm(nsta=2,itemt=c(1,2)) # gaussian and binary items
fit1 <- fitdmm(dat=speed,dmm=mod)

options(digits=3)
cat(fit1$mod$pars,sep=",")

pars <- c(1,0.916,0.084,0.101,0.899,6.39,0.24,0.098,0.902,5.52,0.202,0.472,0.528,1,0)

mod$pars <- pars

options(digits=12)
loglike(speed, mod)
# $logl (for above rounded!!! parameter values) 
# [1] -296.115107102


setwd("/Users/ivisser/Documents/projects/depmixProject/depmixNew/code/depmix/trunk/")

source("depmixS4.R")
source("classes.R")
source("hmModel.R")
source("lystig.R")
source("fb.R")

load("speed.Rda")

# 
# TEST 1: speed data model with optimal parameters, compute the likelihood
# 

pars <- c(1,0.916,0.084,0.101,0.899,6.39,0.24,0.098,0.902,5.52,0.202,0.472,0.528,1,0)

rModels <- list(
  list(
    rModel(formula=rt~1,data=speed,family=gaussian(),pstart=c(5.52,.202)),
    rModel(formula=corr~1,data=speed,family=multinomial(),pstart=c(0.472,0.528))),
  list(
    rModel(formula=rt~1,data=speed,family=gaussian(),pstart=c(6.39,.24)),
    rModel(formula=corr~1,data=speed,family=multinomial(),pstart=c(.098,.902)))
)

trstart=c(0.899,0.101,0.084,0.916)

instart=c(0,1)

mod <- depmix(rModels=rModels,data=speed,trstart=trstart,instart=instart)

ll <- logLik(mod)
ll.fb <- logLik(mod,method="fb")

logl <- -296.115107102 # see above

cat("Test 1: ", all.equal(c(ll),logl,check.att=FALSE), "(loglike of speed data) \n")


# 
# TEST 2
# 
# To check the density function for the multinomial responses with a covariate
# test a model with a single state, which should be identical to a glm
# first fit a model without covariate
# 

invlogit <- function(lp) {exp(lp)/(1+exp(lp))}

acc <- glm(corr~1,data=speed,family=binomial)

p1 <- invlogit(coef(acc)[1])
p0 <- 1-p1

rmod <- list(list(rModel(formula=corr~1,data=speed,family=multinomial(),pstart=c(p0,p1))))
mod <- depmix(rmod,data=speed,trstart=1,instart=c(1))

# check why this does not work with nstates=1!!!!!!!!

# source("depmixS4.R")
# # check why this does not work with nstates=1!!!!!!!!
# 
# source("classes.R")
# source("depmixS4.R")

ll <- logLik(mod)

dev <- -2*ll

cat("Test 2: ", all.equal(c(dev),acc$deviance),"(loglike of 1-comp glm on acc data) \n")


# 
# TEST 3
# 
# now add the covariate and compute the loglikelihood
# 

acc <- glm(corr~Pacc,data=speed,family=binomial)

p1 <- invlogit(coef(acc)[1])
p0 <- 1-p1

rmod <- list(list(rModel(formula=corr~Pacc,data=speed,family=multinomial(),pstart=c(p0,p1,0,coef(acc)[2]))))
mod <- depmix(rmod,data=speed,trstart=1,instart=1)
ll <- logLik(mod)
dev <- -2*ll

cat("Test 3: ", all.equal(c(dev),acc$deviance),"(same but now with covariate) \n")


# 
# TEST 4
# 
# speed accuracy model again, but now with a covariate on the transition matrix
# 

trstart <- c(0.896,0.104,0.084,0.916)

nstates=2
trstart <- matrix(trstart,nstates,byrow=TRUE)

rModels <- list(
  list(
	rModel(formula=rt~1,data=speed,family=gaussian(),pstart=c(5.52,.2)),
	rModel(formula=corr~1,data=speed,family=multinomial())),
  list(
	rModel(formula=rt~1,data=speed,family=gaussian(),pstart=c(6.39,.24)),
	rModel(formula=corr~1,data=speed,family=multinomial(),pstart=c(.098,.902)))
)

trstart=c(0.896,0.104,0.084,0.916)
instart=c(0,1)

trstart=c(trstart[1:2],0,0,trstart[3:4],0,0)
mod <- depmix(rModels=rModels,data=speed,transition=~Pacc,trstart=trstart,instart=instart,ntimes=439)

ll <- logLik(mod)

cat("Test 4a: ", all.equal(c(ll),speedll), "(speed data with cov on transition matrix, value zero) \n")

# start vector with all zeroes, hence prob=F, ie trans probs are specified as mlogit pars and
# not as probabilities
trstart=rep(0,8)
mod1 <- depmix(rModels=rModels,data=speed,transition=~Pacc,instart=instart,trstart=trstart,prob=F) 
ll1 <- logLik(mod1)
# ... which should be identical to the default values chosen if trstart is not provided
mod2 <- depmix(rModels=rModels,data=speed,transition=~Pacc,instart=instart,prob=F) 
ll2 <- logLik(mod2)

cat("Test 4b: ", all.equal(ll1,ll2), "() \n")

rModels <- list(
  list(
	rModel(formula=rt~1,data=speed,family=gaussian(),pstart=c(5.52,.2)),
	rModel(formula=corr~1,data=speed,family=multinomial())),
  list(
	rModel(formula=rt~1,data=speed,family=gaussian(),pstart=c(6.39,.24)),
	rModel(formula=corr~1,data=speed,family=multinomial(),pstart=c(.098,.902)))
)

trstart=c(0.896,0.104,0.084,0.916)
instart=c(0,1)

trstart=c(trstart[1:2],0,0.01,trstart[3:4],0,0.01)

mod <- depmix(rModels=rModels,data=speed,transition=~Pacc,trstart=trstart,instart=instart)

ll <- logLik(mod)

cat("Test 4c: ll is now larger than speedll, ie ll is better due to introduction of a covariate \n")
cat("Test 4c: ", ll,"\t", speedll, "\n")


# 
# test getpars and setpars
# 

rModels <- list(
  list(
	rModel(formula=rt~1,data=speed,family=gaussian(),pstart=c(5.52,.2)),
	rModel(formula=corr~1,data=speed,family=multinomial())),
  list(
	rModel(formula=rt~1,data=speed,family=gaussian(),pstart=c(6.39,.24)),
	rModel(formula=corr~1,data=speed,family=multinomial(),pstart=c(.098,.902)))
)

trstart=c(0.896,0.104,0.084,0.916)
instart=c(0,1)

mod <- depmix(rModels=rModels,data=speed,transition=~1,trstart=trstart,instart=instart)

mod1 <- setpars(mod,getpars(mod))

cat("Test 5: ", all.equal(getpars(mod),getpars(mod1)), "(getpars and setpars) \n")



# 
# SPECIFYING MODELS THE EASY WAY
# 


mod <- depmix(rt~1,data=speed,nstates=2)




