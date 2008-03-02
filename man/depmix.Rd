
\name{depmix}

\alias{depmix}

\title{Dependent Mixture Model Specifiction}

\description{
	
	\item{depmix}{\code{depmix} creates an object of class \code{depmix}, a
	dependent mixture model.}
	
}
\usage{
	depmix(response,data=NULL,nstates,transition,family=gaussian(),prior=~1,initdata=NULL,
		respstart=NULL,trstart=NULL,instart=NULL,ntimes=NULL, ...)
		
}
\arguments{
	\item{response}{The response to be modeled; either a formula or a list 
		of formulae in the multivariate case.}
	
	\item{data}{An optional data.frame to interpret the variables in
		response.}
	
	\item{nstates}{The number of states of the dependent mixture model.}

	\item{transition}{A formula specifying the model for the transition
		matrix; the transition matrix is modeled as a multinomial logistic 
		with possible covariates.}
	
	\item{modname}{A character string with the name of the model, good when 
	fitting many models. Components of mixture models keep their own names. 
	Names are printed in the summary. Boring default names are provided.}
	
	\item{fixed}{A vector of length the number of parameters of the model
	idicating whether parameters are fixed (0) or not (>0).  This may be
	identical to conpat (see below).}
	
	\item{stval}{Start values of the parameters.  These will be random if
	not specified.  Start values must be specified (for all parameters) if
	there are fixed parameters.}
	
	\item{conrows}{Argument \code{conrows} can be used to specify general
	constraints between parameters. See details below.}
	
	\item{conpat}{Argument \code{conpat} can be used to specify fixed
	parameters and equality constraints.  It can not be used in conjuction
	with fixed.  See details below.}
	
	\item{tdfix,tdst}{The first is a logical vector indicating (with 1's) 
	which parameters are dependent on covariates (it should have length npars). 
	tdst provides the starting values for the regression parameters. 
	Using tdcov=TRUE in fitdepmix will actually fit the regression parameters. 
	The covariate itself has to be specified in the data as 
	"covariate" (see help on markovdata) and should be scaled to 0-1.}

	\item{linmat}{A complete matrix of linear constraints.  This argument
	is intended for internal use only, it is used by the fit routine to
	re-create the model with the fitted parameter values.  Warning: use of
	this argument results in complete replacement of the otherwise created
	matrix A, which contains e.g. sum contraints for transition matrix
	parameters.  If \code{linmat} is provided, make sure it is correct,
	otherwise strange results may occur in fitting models.}
	
	\item{snames}{Names for the states may be provided in snames. 
	Defaults are State1, State2 etc. These are printed in the summary.}

	\item{inames}{Names for items may be provided in inames. 
	Defaults are Item1, Item2 etc. They are printed in the summary.}
		
	\item{depmix}{Object of class \code{depmix}.}
		
	\item{precision}{Precision sets the number of digits to be printed in the
	summary functions.}
	
	\item{se}{Vector with standard errors, these are passed on from the 
	summary.fit function if and when ses are available.}
	
	\item{specs,...}{Internal use.}
	
	\item{object}{An object of class \code{depmix}.}
}

\details{ 

	The function \code{depmix} creates an object of class \code{depmix} and sets
	random initial parameter values if these are not provided.  Even though
	\code{depmix} is not a mixture of Markov models, the mixture parameter is
	is included in the parameter vector.  This is important when specifying
	constraints.  Parameters are ordered as follows: the first parameter(s)
	are the mixing proportions of the mixture of Markov and/or latent class
	models.  I.e., when a single latent class model or a single Markov
	chain is fitted, this mixture proportion has value 1.0 and is it is
	fixed in estimation.  After the mixing proportions, the next parameters
	in the parameter vector are the transition matrix parameters, the
	square of nstates in row-major order.  That is, first the transition
	probabilities from state 1 to all the other states are given, then the
	probabilities from state 2 to all the other states etc.  Next are the
	observation matrix parameters.  These are provided consecutively for
	each state/class.  Ie a trichtomous item model with two states has 6
	observation parameters; the first three are the probabilities of
	observing category 1, 2 and 3 respectively in state 1 (which sum to
	one), and then similarly for state 2.  As another example: suppose we
	have model for one binary item and one gaussian item, in that order, we
	would have 4 observation parameters for each state, first the
	probabilities of observing a symbol from category 1 or 2 in state 1,
	the two parameters, the mean and standard deviation for state 1, and
	then the same state 2 (see the example in fitdepmix with data from rudy).
	Finally the initial state probabilities are provided, in the order of
	the states.  In the case of a latent class model or a finite mixture
	model, these parameters are usually denote as the mixture proportions.
	
	Linear constraints can be set using arguments \code{conrows} and
	\code{conpat}.  \code{conrows} must be contain nc by npars values, in
	row major order, with nc the number of contraints to be specified.
	\code{conrows} is used to define general linear constraints.  A row of
	\code{conrows} must contain the partial derivatives of a general linear
	constraint with respect to each of the parameters.  Suppose we want the
	constraint x1 -2*x2=0, one row of conrows should contain a 1 in
	position one and -2 in position and zeroes in the remaining positions.
	In the function \code{mixdepmix} \code{conrows} is understood to specify
	linear constraints on the mixing proportions only.  As a consequence,
	it is not possible to easily constrain parameters between components of
	a mixture model.

	\code{conpat} can be used as a shortcut for both fixed and conrows.  It
	must be a single vector of length npars contaning 0's (zeroes) for
	fixed parameters, 1's (ones) for free parameters and higher numbers for
	possibly equality constrained parameters.  E.g.
	\code{conpat=c(1,1,0,2,2,3,3,3)} would indicate that pars 1 and 2 are
	freely estimated, par 3 is fixed at its startvalue (which must be
	provided in this case), par 4 and 5 are to estimated equal and pars 6,
	7 and 8 are also to be estimated equal.

}
\value{
	
\code{depmix} returns an object of class \code{depmix} which has its own summary
method.  This will print the parameter values, itemtypes, number of (free)
parameters, and the number of states.  There is no print method.  Using
print will print all fields of the model which is a list of the following:

	\item{modname}{See above.}
	
	\item{nstates}{See above}

	\item{snames}{See above.}

	\item{nitems}{The number of items(=length(itemtypes)).}
	
	\item{itemtypes}{See above.}
	
	\item{inames}{See above.}

	\item{npars}{The total parameter count of the model.}
	
	\item{nparstotal}{The total number of parameters of when the covariate 
	parameters are included.}
	
	\item{freepars}{The number of freely estimated parameters (it is
	computed as sum(as.logical(fixed))-rank(qr(A)).}
	
	\item{freeparsnotd}{The number of freely estimated parameters (it is
	computed as sum(as.logical(fixed))-rank(qr(A)); this version without
	the covariate parameters.}

	\item{pars}{A vector of length npars containing parameter values.}
	
	\item{fixed}{\code{fixed} is a (logical) vector of length npars
	specifying which parameters are fixed and which are not.}
	
	\item{A}{The matrix A contains the general linear constraints of the
	model.  nrow(A) is the number of linear constraints.  A starts with a
	number of rows for the sum constraints for the transition, observation
	and initial state parameters, after which the user provided constraints
	are added.}

	\item{bu,bl}{bu and bl represent the upper and lower bounds of the
	parameters and the constraints.  These vectors are each of length npars
	+ nrow(A).}
	
	\item{bllin,bulin}{The lower and upper bounds of the linear
	constraints.}

	\item{td,tdin,tdtr,tdob,tdfit}{Logicals indicating whehter there
	covariates, in which parameters they are, and whether they are
	estimated or not (the latter is used to decide whether to print those
	values or not).}

	\item{st}{Logical indicating whether the model has user specified
	starting values.}
	
\code{lca} returns an object of class \code{depmix, lca}, and is otherwise
identical to a \code{depmix} object. The only difference is that the
transition matrix parameters are irrelevant, and
consequently they are not printed in the summary function. 

}

\author{Ingmar Visser \email{i.visser@uva.nl}} 

\seealso{\code{\link[depmix]{mixdmm}} on defining mixtures of \code{dmm}'s,
\code{\link[depmix]{mgdmm}} for defining multi group models, and
\code{\link[depmix]{generate}} for generating data from models. 
}

\references{
	
	On hidden Markov models: Lawrence R. Rabiner (1989).  A tutorial on
	hidden Markov models and selected applications in speech recognition.
	\emph{Proceedings of IEEE}, 77-2, p.  267-295.
	
	On latent class models: A. L. McCutcheon (1987).  \emph{Latent class
	analysis}.  Sage Publications.
	
}

\examples{

# create a 2 state model with one continuous and one binary response
# with start values provided in st
st <- c(1,0.9,0.1,0.2,0.8,2,1,0.7,0.3,5,2,0.2,0.8,0.5,0.5)
mod <- dmm(nsta=2,itemt=c(1,2), stval=st)
summary(mod)

# 2 class latent class model with equal conditional probabilities in each class
stv=c(1,rep(c(0.9,0.1),5),rep(c(0.1,0.9),5),0.5,0.5)
# here the conditional probs of the first item are set equal to those in
# the subsequent items
conpat=c(1,rep(c(2,3),5),rep(c(4,5),5),1,1)
lc=lca(ncl=2,itemtypes=rep(2,5),conpat=conpat,stv=stv)
summary(lc)

}

\keyword{models}


