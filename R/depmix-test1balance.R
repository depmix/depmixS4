
# 
# Started by Ingmar Visser 16-10-2007
# 
# Usage: go to trunk directory and source this file in R, if the program
# still works it should return TRUE at every test (or make immediate sense
# otherwise)

# Changes: 

# 
# TEST ntimes argument using balance scale data, ie lca data
# 

# old depmix results

setwd("/Users/ivisser/Documents/projects/depmixProject/depmixNew/code/depmix/trunk/")

.libPaths(new="/Users/ivisser/Library/R/library/")
library(depmix)

# balance scale data: distance items only
# nit=5
# nt=rep(1,472)
# itt=rep("cat",5)
# dat=scan("afstand")
# dat[which(dat==2)]=0
# balance <- markovdata(dat,item=itt,ntim=nt,dname="Balance scale data",
# 	ina=c("d1","d2","d3","d4","d5"))
# 
# save(balance,file="balance.rda",ascii=TRUE)

load("balance.rda")

summary(balance)

bmod <- lca(nc=2,itemt=rep(2,5))
fit1 <- fitdmm(balance,bmod,meth="npsol")

pars <- c(1,
	0.886,0.114,
	0.967,0.0326,
	0.877,0.123,
	0.9,0.1,
	0.934,0.0657,
	0.0903,0.91,
	0.0814,0.919,
	0.0559,0.944,
	0.0795,0.92,
	0.0673,0.933,
	0.262,0.738)

bmod <- lca(nc=2,itemt=rep(2,5),st=pars)

options(digits=12)
loglike(balance,bmod)

# $logl
# [1] -898.39684242

# compare with optimized loglike (non-rounded par values)
# $logl
# [1] -898.396685164


# load depmixNew

setwd("/Users/ivisser/Documents/projects/depmixProject/depmixNew/code/depmix/trunk/")

source("depmixS4.r")
source("classes.r")
source("hmModel.R")
source("fithmModel.R")
source("llratio.R")
source("lystig.R")
source("fb.R")


# change data to data.frame format
load("balance.rda")
cnames <- colnames(balance)
balance <- matrix(balance,ncol=5)
hist(rowSums(balance))

balance <- data.frame(balance)
colnames(balance) <- cnames

hist(rowSums(balance))

pars <- c(1,
	0.886,0.114,
	0.967,0.0326,
	0.877,0.123,
	0.9,0.1,
	0.934,0.0657,
	0.0903,0.91,
	0.0814,0.919,
	0.0559,0.944,
	0.0795,0.92,
	0.0673,0.933,
	0.262,0.738)

# above par values should result in this loglikelihood
# $logl
# [1] -898.39684242

# 
# define latent class model
# 


rModels <- list(
  list(
	rModel(formula=d1~1,data=balance,family=multinomial(),pstart=c(0.9,0.1)),
	rModel(formula=d2~1,data=balance,family=multinomial(),pstart=c(0.9,0.1)),
	rModel(formula=d3~1,data=balance,family=multinomial(),pstart=c(0.9,0.1)),
	rModel(formula=d4~1,data=balance,family=multinomial(),pstart=c(0.9,0.1)),
	rModel(formula=d5~1,data=balance,family=multinomial(),pstart=c(0.9,0.1))),
  list(
	rModel(formula=d1~1,data=balance,family=multinomial(),pstart=c(0.1,0.9)),
	rModel(formula=d2~1,data=balance,family=multinomial(),pstart=c(0.1,0.9)),
	rModel(formula=d3~1,data=balance,family=multinomial(),pstart=c(0.1,0.9)),
	rModel(formula=d4~1,data=balance,family=multinomial(),pstart=c(0.1,0.9)),
	rModel(formula=d5~1,data=balance,family=multinomial(),pstart=c(0.1,0.9)))
)

trstart=c(1,0,0,1)

instart=c(0.262,0.738)



# ntimes is added as an argument as the attribute ntimes for speed is different from this
mod <- depmix(rModels=rModels,data=balance,trstart=trstart,instart=instart,ntimes=rep(1,472))

logLik(mod)

source("fithmModel.R")

pars <- getpars(mod)
fixed <- c(1,0,1,1,1,1,rep(c(1,0),10))

mod1 <- fit(mod,fixed=fixed)

logLik(mod1)


# 
# example with a covariate on the initial state probs
# 

sumscores <- rowSums(balance)

set.seed(12)

age <- (sumscores+12)/2+rnorm(472,sd=1)
cor(age,sumscores)

balance$age <- age

instart=c(0.262,0.738,0,0.01)

# ntimes is added as an argument as the attribute ntimes for speed is different from this
mod <- depmix(rModels=rModels,init=~age,instart=instart,data=balance,
	trstart=trstart,ntimes=rep(1,472))

logLik(mod)

# logLik is now better than the original -898.3968 due to introduction of
# the covariate on the init probs











