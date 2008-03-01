
# 
# Started by Ingmar Visser 16-10-2007
# 
# Usage: go to trunk directory and source this file in R, if the program
# still works it should return TRUE at every test (or make immediate sense
# otherwise)

# Changes: 

# 
# TEST ntimes argument using discrimination data
# 

# old depmix results

.libPaths(new="/Users/ivisser/Library/R/library/")
library(depmix)

setwd("/Users/ivisser/Documents/projects/depmixProject/depmixNew/code/depmix/trunk/")

# all together
ntimes=c(11,17,18,11,25,15,14,10,13,32,10,10,11,10,12,15,14,10,15,23,10,
13,15,17,9,11,26,20,18,13,9,10,11,9,9,18,20,22,45,13,25,11,
11,15,14,12,23,11,47,27,16,10,12,15,10,14,18,17,10,16,11,23,11,
18,12,10,11,16,10,9,11,10,10,9,13,11,10,10,11,11,10,10,37,9,
11,10,10,10,9,10,10,10,10,10,16,13,10,12,13,9,11,9,15,11,10,
13,12,9,12,14,22,14,15,10,10,13,10,9,12,11,9,10,10,11,10,10,
12,9,9,16,12,10,11,11,12,10,9,14,10,11,10,10,9,10,10,10,10,
10,14,11,12,10,10,9,11,10,14,25,23,39,29,29,26,23,21,29,28,24,
33,27,23,41,20,22,31,29,48,22,44,42,19,30,35,39,40,43,38,17,27,
22,47,42)
dat=scan("discrimination")

dat=replace(dat,which(dat==2),0)

discrimination=markovdata(dat=dat,item="cat",ntimes=ntimes,dname="discrimination",inames="correct")


# all or none model with error prob in the learned state
fixed = c(0,0,0,1,1,1,1,0,0,0,0)
stv = c(1,1,0,0.03,0.97,0.1,0.9,0.5,0.5,0,1)
allor <- dmm(nstates=2,itemtypes=2,fixed=fixed,stval=stv,modname="All-or-none")

fit1 <- fitdmm(dat=discrimination,dmm=allor,meth="npsol")

fmod <- fit1$mod
options(digits=4)
cat(fmod$pars)

pars=c(1,1,0,0.118,0.882,0.0501,0.9499,0.5,0.5,0,1)

fmod$pars <- pars

options(digits=12)
loglike(discrimination,fmod)

# above (rounded!) par values result in this loglike
# $logl
# [1] -1670.77403661

# 
# new depmix version results
# 

setwd("/Users/ivisser/Documents/projects/depmixProject/depmixNew/code/depmix/trunk/")

source("depmixS4.R")
source("classes.R")
source("hmModel.R")


# change data to data.frame format

ntimes=c(11,17,18,11,25,15,14,10,13,32,10,10,11,10,12,15,14,10,15,23,10,
13,15,17,9,11,26,20,18,13,9,10,11,9,9,18,20,22,45,13,25,11,
11,15,14,12,23,11,47,27,16,10,12,15,10,14,18,17,10,16,11,23,11,
18,12,10,11,16,10,9,11,10,10,9,13,11,10,10,11,11,10,10,37,9,
11,10,10,10,9,10,10,10,10,10,16,13,10,12,13,9,11,9,15,11,10,
13,12,9,12,14,22,14,15,10,10,13,10,9,12,11,9,10,10,11,10,10,
12,9,9,16,12,10,11,11,12,10,9,14,10,11,10,10,9,10,10,10,10,
10,14,11,12,10,10,9,11,10,14,25,23,39,29,29,26,23,21,29,28,24,
33,27,23,41,20,22,31,29,48,22,44,42,19,30,35,39,40,43,38,17,27,
22,47,42)
dat=scan("discrimination")
dat=replace(dat,which(dat==2),0)
disc <- data.frame(dat)

colnames(disc) <- "acc"

pars=c(1,1,0,0.118,0.882,0.0501,0.9499,0.5,0.5,0,1)

# 
# define all-or-none model
# 

source("depmixS4.R")

rModels <- list(
  list(
	rModel(formula=acc~1,data=disc,family=multinomial(),pstart=c(0.0501,0.9499))),
  list(
	rModel(formula=acc~1,data=disc,family=multinomial(),pstart=c(0.5,0.5)))
)

trstart=c(1,0,0.118,0.882)

instart=c(0,1)

# ntimes is added as an argument as the attribute ntimes for speed is different from this
mod <- depmix(rModels=rModels,data=disc,trstart=trstart,instart=instart,ntimes=ntimes)

mod

source("depmixS4.R")

options(digits=12)
logLik(mod,meth="fb")

# result
# [1] -1670.77403661

# 
# TEST SUCCESFULL!!!!!!!
# 




