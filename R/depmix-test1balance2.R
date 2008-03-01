
# 
# Started by Ingmar Visser 26-2-2008
# 
# Usage: go to trunk directory and source this file in R, if the program
# still works it should return TRUE at every test (or make immediate sense
# otherwise)

# Changes: 

# 
# BALANCE SCALE data example with age as covariate on class membership
# 

# old depmix results

.libPaths(new="/Users/ivisser/Library/R/library/")
library(depmix)

setwd("/Users/ivisser/Documents/projects/depmixProject/depmixNew/code/depmix/trunk/")

balance=read.table("balance.txt",head=T)

# the data set consists of four distance items on a balance scale task
# the data frame contains 11 columns with covariates sex age in days and age in years,
# and two versions of the four items; t1-t4 are trichotomous items with 1, 2, and 3 
# meaning left, balance, and right; d1-d4 are dichotomously scored items with 0=incorrect
# and 1=correct

# to get to grips with the data, first study the histogram of the sumscores of the
# dichotomous items

sums <- apply(balance[,8:11],1,sum)
hist(sums)


source("depmixS4.r")
source("classes.r")
source("hmModel.R")
source("fithmModel.R")
source("llratio.R")
source("lystig.R")
source("fb.R")


# now fit some latent class models

rMods <- list(
  list(
	rModel(formula=d1~1,data=balance,family=multinomial()),
	rModel(formula=d2~1,data=balance,family=multinomial()),
	rModel(formula=d3~1,data=balance,family=multinomial()),
	rModel(formula=d4~1,data=balance,family=multinomial())),
  list(
	rModel(formula=d1~1,data=balance,family=multinomial()),
	rModel(formula=d2~1,data=balance,family=multinomial()),
	rModel(formula=d3~1,data=balance,family=multinomial()),
	rModel(formula=d4~1,data=balance,family=multinomial()))
)

trstart=c(1,0,0,1)
instart=c(0.5,0.5)

# ntimes is added as an argument as the attribute ntimes for speed is different from this
mod <- depmix(rModels=rMods,data=balance,trstart=trstart,instart=instart,ntimes=rep(1,nrow(balance)))

pars <- getpars(mod)
fixed <- c(1,0,1,1,1,1,rep(c(1,0),8))

mod1 <- fit(mod,fixed=fixed)

logLik(mod1)

AIC(mod1)
BIC(mod1)


logit <- function(p) {
	log(p/(1-p))
}

invlogit <- function(x) {
	exp(x)/(1+exp(x))
}




# trichotome data
dat3 <- dat[,c(4:7)]
dat3 <- markovdata(dat3,nt=rep(1,nrow(dat3)),itemt=c(3,3,3,3))


# 
# TRICHOTOME DATA MODELLEN
# 

lc1 <- lca(nc=1,itemt=c(3,3,3,3))
fit1 <- fitdmm(dat3,lc1)

lc2 <- lca(nc=2,itemt=c(3,3,3,3))
fit2 <- fitdmm(dat3,lc2)

lc3 <- lca(nc=3,itemt=c(3,3,3,3))
fit3.2 <- fitdmm(dat3,lc3)

lc4 <- lca(nc=4,itemt=c(3,3,3,3))
fit4 <- fitdmm(dat3,lc4)



Model:  1 class model  fitted at  Mon Feb 18 11:48:00 2008  
Optimization information, method is  donlp 
 Iterations:  37
 Inform:  KT-conditions satisfied, no further correction computed  (look up the respective manuals for more information.)

 Loglikelihood of fitted model:  -2227.966
 AIC:  4471.932
 BIC:  4509.196
 Number of observations (used in BIC):  779
 Fitted model 
 Model:  1 class model  
 Number of parameters:  15  
 Free parameters:       8  
 Number of states:      1 
 Number of items:       4 
 Item types:            3 3 3 3 

 Parameter values, observation parameters 

	   Item1,p 1 Item1,p 2 Item1,p 3 Item2,p 1 Item2,p 2 Item2,p 3 Item3,p 1 Item3,p 2
Class1     0.659     0.322     0.019     0.018     0.302     0.680     0.653     0.309
se         0.017     0.017     0.005     0.005     0.016     0.017     0.017     0.017
t         38.760    19.244     3.911     3.776    18.344    40.720    38.322    18.680
	   Item3,p 3 Item4,p 1 Item4,p 2 Item4,p 3
Class1     0.037     0.017     0.294     0.689
se         0.007     0.005     0.016     0.017
t          5.488     3.636    18.010    41.577


 Parameter values, unconditional (class) probabilities 

	Class1
val      1
se       0
t       NA


Model:  2 class model  fitted at  Mon Feb 18 11:50:22 2008  
Optimization information, method is  donlp 
 Iterations:  30
 Inform:  KT-conditions satisfied, no further correction computed  (look up the respective manuals for more information.)

 Loglikelihood of fitted model:  -1326.523
 AIC:  2687.046
 BIC:  2766.232
 Number of observations (used in BIC):  779
 Fitted model 
 Model:  2 class model  
 Number of parameters:  31  
 Free parameters:       17  
 Number of states:      2 
 Number of items:       4 
 Item types:            3 3 3 3 

 Parameter values, observation parameters 

	   Item1,p 1 Item1,p 2 Item1,p 3 Item2,p 1 Item2,p 2 Item2,p 3 Item3,p 1 Item3,p 2
Class1     0.036     0.919     0.045     0.037     0.892     0.070     0.030     0.895
se         0.015     0.020     0.014     0.012     0.023     0.021     0.013     0.021
t          2.393    46.329     3.293     3.021    38.103     3.392     2.270    42.384
Class2     0.938     0.055     0.008     0.009     0.037     0.954     0.933     0.047
se         0.012     0.011     0.004     0.004     0.009     0.010     0.013     0.011
t         78.442     4.879     1.875     2.185     4.086    95.947    73.884     4.270
	   Item3,p 3 Item4,p 1 Item4,p 2 Item4,p 3
Class1     0.075     0.041     0.893     0.067
se         0.017     0.013     0.024     0.022
t          4.354     3.151    36.567     3.083
Class2     0.020     0.006     0.026     0.968
se         0.006     0.003     0.007     0.008
t          3.215     1.718     3.423   117.896


 Parameter values, unconditional (class) probabilities 

	Class1 Class2
val  0.309  0.691
se   0.018  0.018
t   17.535 39.123

Model:  3 class model  fitted at  Mon Feb 18 11:54:10 2008  
Optimization information, method is  donlp 
 Iterations:  104
 Inform:  KT-conditions satisfied, no further correction computed  (look up the respective manuals for more information.)

 Loglikelihood of fitted model:  -1243.951
 AIC:  2539.902
 BIC:  2661.01
 Number of observations (used in BIC):  779
 Fitted model 
 Model:  3 class model  
 Number of parameters:  49  
 Free parameters:       26  
 Number of states:      3 
 Number of items:       4 
 Item types:            3 3 3 3 

 Parameter values, observation parameters 

	   Item1,p 1 Item1,p 2 Item1,p 3 Item2,p 1 Item2,p 2 Item2,p 3 Item3,p 1 Item3,p 2
Class1     0.938     0.056     0.006     0.009     0.038     0.953     0.935     0.048
se         0.012     0.011     0.004     0.004     0.009     0.010     0.012     0.011
t         79.235     4.973     1.536     2.187     4.133    96.004    74.929     4.335
Class2     0.075     0.105     0.819     0.676     0.000     0.324     0.087     0.090
se         0.084       NaN     0.072     0.317     0.387     0.180     0.084     0.103
t          0.893        NA    11.377     2.133     0.000     1.802     1.037     0.877
Class3     0.036     0.960     0.004     0.000     0.939     0.061     0.026     0.935
se         0.015     0.016     0.005       NaN     0.007     0.020     0.013     0.017
t          2.377    59.545     0.775        NA   133.115     3.016     2.057    55.569
	   Item3,p 3 Item4,p 1 Item4,p 2 Item4,p 3
Class1     0.017     0.004     0.026     0.970
se         0.006     0.003     0.007     0.008
t          2.873     1.393     3.458   121.763
Class2     0.823     0.733     0.173     0.094
se         0.123     0.152     0.054     0.142
t          6.719     4.833     3.198     0.667
Class3     0.039     0.004     0.930     0.065
se         0.012     0.004     0.022     0.022
t          3.337     1.001    41.420     2.959


 Parameter values, unconditional (class) probabilities 

	Class1 Class2 Class3
val  0.689  0.017  0.294
se   0.018  0.005  0.017
t   39.127  3.331 16.941

Model:  4 class model  fitted at  Mon Feb 18 11:59:58 2008  
Optimization information, method is  donlp 
 Iterations:  80
 Inform:  KT-conditions (relaxed) satisfied, singular point  (look up the respective manuals for more information.)

 Loglikelihood of fitted model:  -1220.541
 AIC:  2511.081
 BIC:  2674.112
 Number of observations (used in BIC):  779
 Fitted model 
 Model:  4 class model  
 Number of parameters:  69  
 Free parameters:       35  
 Number of states:      4 
 Number of items:       4 
 Item types:            3 3 3 3 

 Parameter values, observation parameters 

	   Item1,p 1 Item1,p 2 Item1,p 3 Item2,p 1 Item2,p 2 Item2,p 3 Item3,p 1 Item3,p 2
Class1     0.978     0.020     0.002     0.007     0.018     0.975     0.981     0.011
se         0.032     0.027     0.006     0.005     0.014     0.016     0.037     0.028
t         30.972     0.765     0.338     1.325     1.328    59.916    26.512     0.410
Class2     0.022     0.978     0.000     0.000     0.972     0.028     0.017     0.953
se         0.013     0.072     0.073       NaN     0.018     0.020     0.011     0.018
t          1.651    13.651     0.000        NA    54.756     1.437     1.493    52.043
Class3     0.087     0.079     0.834     0.744     0.000     0.256     0.087     0.085
se         0.089       NaN     0.043     0.245     0.261     0.131     0.090     0.084
t          0.980        NA    19.473     3.038     0.000     1.953     0.961     1.012
Class4     0.499     0.448     0.053     0.025     0.283     0.692     0.440     0.441
se         0.295     0.307     0.040     0.028     0.248     0.249     0.315     0.294
t          1.690     1.458     1.335     0.878     1.143     2.782     1.396     1.502
	   Item3,p 3 Item4,p 1 Item4,p 2 Item4,p 3
Class1     0.007     0.004     0.017     0.979
se         0.011     0.003     0.010     0.010
t          0.692     1.412     1.695    98.374
Class2     0.030     0.005     0.977     0.018
se         0.013     0.005     0.028     0.028
t          2.359     1.022    34.440     0.641
Class3     0.828     0.836     0.164     0.000
se         0.118     0.268     0.060     0.273
t          7.012     3.126     2.707     0.000
Class4     0.119     0.000     0.177     0.823
se         0.050       NaN     0.182     0.194
t          2.354        NA     0.973     4.250


 Parameter values, unconditional (class) probabilities 

	Class1 Class2 Class3 Class4
val  0.616  0.270  0.015  0.099
se   0.088  0.023  0.004  0.074
t    7.001 11.687  3.667  1.331



# 
# 3-class model is best according to BIC
# 

# get the posteriors from it

dat$post <- as.factor(fit3$post$states[[1]][,2])

levels(dat$post) <- c("rule2","trans","rule1")

> lm(ageyears~post,dat=dat)

Call:
lm(formula = ageyears ~ post, data = dat)

Coefficients:
(Intercept)    posttrans    postrule1  
	12.7908      -0.3293      -3.6279  

> plot(ageyears~post,dat=dat)
> summary(multinom(post~ageyears,dat=dat))
# weights:  9 (4 variable)
initial  value 855.818973 
iter  10 value 398.743386
final  value 393.199188 
converged
Call:
multinom(formula = post ~ ageyears, data = dat)

Coefficients:
	  (Intercept)    ageyears
trans   -3.190456 -0.04319348
rule1    5.233625 -0.56814085

Std. Errors:
	  (Intercept)   ageyears
trans   1.2964214 0.10150279
rule1   0.4683135 0.04524304

Residual Deviance: 786.3984 
AIC: 794.3984 






# 
# DICHOTOME DATA MODELLEN
# 


setwd("/Users/ivisser/Documents/projects/depmixProject/lcavoorbeeld/")

dat=read.table("lca.txt",head=T)

# trichotome data
dat <- dat[,c("d2_di","d3_bi","d4_di","d5_di","ageyears")]

Model:  1 class model  fitted at  Fri Jan 25 14:48:04 2008  
Optimization information, method is  donlp 
 Iterations:  10
 Inform:  KT-conditions satisfied, no further correction computed  (look up the respective manuals for more information.)

 Loglikelihood of fitted model:  -1973.626
 AIC:  3955.253
 BIC:  3973.885
 Number of observations (used in BIC):  779
 Fitted model 
 Model:  1 class model  
 Number of parameters:  11  
 Free parameters:       4  
 Number of states:      1 
 Number of items:       4 
 Item types:            2 2 2 2 

 Parameter values, observation parameters 

	   Item1,p 1 Item1,p 2 Item2,p 1 Item2,p 2 Item3,p 1 Item3,p 2 Item4,p 1 Item4,p 2
Class1     0.341     0.659     0.320     0.680     0.347     0.653     0.311     0.689
se         0.017     0.017     0.017     0.017     0.017     0.017     0.017     0.017
t         20.098    38.761    19.131    40.720    20.328    38.322    18.737    41.576


 Parameter values, unconditional (class) probabilities 

	Class1
val      1
se       0
t       NA




Model:  2 class model  fitted at  Fri Jan 25 14:47:02 2008  
Optimization information, method is  donlp 
 Iterations:  18
 Inform:  KT-conditions satisfied, no further correction computed  (look up the respective manuals for more information.)

 Loglikelihood of fitted model:  -1083.036
 AIC:  2184.073
 BIC:  2225.995
 Number of observations (used in BIC):  779
 Fitted model 
 Model:  2 class model  
 Number of parameters:  23  
 Free parameters:       9  
 Number of states:      2 
 Number of items:       4 
 Item types:            2 2 2 2 

 Parameter values, observation parameters 

	   Item1,p 1 Item1,p 2 Item2,p 1 Item2,p 2 Item3,p 1 Item3,p 2 Item4,p 1 Item4,p 2
Class1     0.965     0.035     0.929     0.071     0.970     0.030     0.934     0.066
se         0.015     0.015     0.021     0.021     0.013     0.013     0.022     0.022
t         64.440     2.314    43.618     3.350    74.558     2.270    42.164     2.982
Class2     0.062     0.938     0.047     0.953     0.067     0.933     0.031     0.969
se         0.012     0.012     0.010     0.010     0.013     0.013     0.008     0.008
t          5.123    77.910     4.680    95.892     5.220    72.868     3.814   118.626


 Parameter values, unconditional (class) probabilities 

	Class1 Class2
val  0.310  0.690
se   0.018  0.018
t   17.449 38.911


Model:  3 class model  fitted at  Fri Jan 25 14:49:39 2008  
Optimization information, method is  donlp 
 Iterations:  27
 Inform:  KT-conditions (relaxed) satisfied, singular point  (look up the respective manuals for more information.)

 Loglikelihood of fitted model:  -1062.374
 AIC:  2152.747
 BIC:  2217.96
 Number of observations (used in BIC):  779
 Fitted model 
 Model:  3 class model  
 Number of parameters:  37  
 Free parameters:       14  
 Number of states:      3 
 Number of items:       4 
 Item types:            2 2 2 2 

 Parameter values, observation parameters 

	   Item1,p 1 Item1,p 2 Item2,p 1 Item2,p 2 Item3,p 1 Item3,p 2 Item4,p 1 Item4,p 2
Class1     0.982     0.018     0.962     0.038     0.981     0.019     0.982     0.018
se         0.038     0.038     0.071     0.071     0.027     0.027     0.137     0.137
t         26.003     0.471    13.548     0.529    36.056     0.699     7.180     0.128
Class2     0.023     0.977     0.024     0.976     0.017     0.983     0.019     0.981
se         0.109     0.109     0.055     0.055     0.159     0.159     0.030     0.030
t          0.211     8.943     0.438    17.879     0.109     6.168     0.647    32.936
Class3     0.473     0.527     0.314     0.686     0.560     0.440     0.201     0.799
se         1.376     1.376     1.155     1.155     1.365     1.365     0.854     0.854
t          0.344     0.383     0.272     0.594     0.410     0.322     0.236     0.935


 Parameter values, unconditional (class) probabilities 

	Class1 Class2 Class3
val  0.283  0.613  0.105
se   0.099  0.381  0.284
t    2.864  1.609  0.369

# posterior classes related to age in years
# weights:  9 (4 variable)
initial  value 855.818973 
iter  10 value 526.053934
final  value 526.041488 
converged
Call:
multinom(formula = post ~ ageyears, data = dat)

Coefficients:
  (Intercept)  ageyears
2   -5.320429 0.5578123
3   -4.356987 0.3028828

Std. Errors:
  (Intercept)   ageyears
2   0.4641045 0.04387059
3   0.6461756 0.05985131

Residual Deviance: 1052.083 
AIC: 1060.083 



Model:  4 class model  fitted at  Fri Jan 25 14:53:45 2008  
Optimization information, method is  donlp 
 Iterations:  38
 Inform:  KT-conditions (relaxed) satisfied, singular point  (look up the respective manuals for more information.)

 Loglikelihood of fitted model:  -1058.019
 AIC:  2154.038
 BIC:  2242.54
 Number of observations (used in BIC):  779
 Fitted model 
 Model:  4 class model  
 Number of parameters:  53  
 Free parameters:       19  
 Number of states:      4 
 Number of items:       4 
 Item types:            2 2 2 2 

 Parameter values, observation parameters 

	   Item1,p 1 Item1,p 2 Item2,p 1 Item2,p 2 Item3,p 1 Item3,p 2 Item4,p 1 Item4,p 2
Class1     0.995     0.005     0.987     0.013     1.000     0.000     0.978     0.022
se         0.042     0.042       NaN       NaN       NaN       NaN     0.078     0.078
t         23.685     0.121        NA        NA        NA        NA    12.507     0.278
Class2     0.700     0.300     0.000     1.000     1.000     0.000     0.388     0.612
se         0.063     0.063     2.732     2.732       NaN       NaN       NaN       NaN
t         11.158     4.786     0.000     0.366        NA        NA        NA        NA
Class3     0.567     0.433     1.000     0.000     0.511     0.489     0.527     0.473
se           NaN       NaN     0.029     0.029     1.122     1.122       NaN       NaN
t             NA        NA    34.293     0.000     0.456     0.436        NA        NA
Class4     0.042     0.958     0.025     0.975     0.033     0.967     0.023     0.977
se         0.009     0.009     0.009     0.009     0.011     0.011     0.007     0.007
t          4.632   106.510     2.722   105.363     2.894    85.135     3.372   143.828


 Parameter values, unconditional (class) probabilities 

	Class1 Class2 Class3 Class4
val  0.263  0.040  0.044  0.654
se   0.016  0.109  0.101  0.016
t   16.785  0.365  0.433 41.634






