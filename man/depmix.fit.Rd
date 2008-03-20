\name{fit}

\docType{method}

\alias{fit}
\alias{depmix.fit}
\alias{fit,depmix-method}

\alias{show,depmix.fitted-method}
\alias{summary,depmix.fitted-method}

\title{Fit depmix models}

\description{

	\code{fit} optimizes parameters of \code{\link{depmix}} models, optionally
	subject to general linear (in)equality constraints.}

\usage{
	
	\S4method{fit}{depmix}(object, fixed=NULL, equal=NULL, conrows=NULL,
		conrows.upper=0, conrows.lower=0, method=NULL,...)
	
	\S4method{summary}{depmix.fitted}(object)
	
}

\arguments{

	\item{object}{An object of class \code{depmix}.}
	
	\item{fixed}{Vector of mode logical indicating which parameters should 
		be fixed.}

	\item{equal}{Vector indicating equality constraints; see details.}

	\item{conrows}{Rows of a general linear constraint matrix; see details.}

	\item{conrows.upper, conrows.lower}{Upper and lower bounds for the
		linear constraints; see details.}

	\item{method}{The optimization method; mostly determined by
		constraints.}

	\item{...}{Further arguments passed on to the optimization methods.}

}

\details{
	The method fits depmix models by the EM algorithm if there are no
	linear constraints on the parameters and if the transition model has
	no covariates. Otherwise the general optimizer donlp is used which 
	handles general linear (in-)equality constraints. 
	
	Three types of constraints can be specified on the parameters: fixed,
	equality, and general constraints. Constraint vectors should be of
	length npar(object). See help on getpars and setpars about the
	ordering of parameters. 
	
	The \code{equal} argument is used to specify equality constraints:
	parameters that get the same integer number in this vector are
	estimated to be equal. Any integers can be used in the vector except 0
	and 1, which indicate fixed and free parameters respectively. 

	Using the donlp optimizer a Newton-Raphson scheme is employed to
	estimate parameters subject to linear constraints by imposing: 
	
			bl <= A*x <= bu,
	
	where x is the parameter vector, bl is a vector of lower bounds, bu
	is a vector of upper bounds, and A is the constraint matrix. 

	The \code{conrows} argument is used to specify rows of A directly, and 
	the conrows.lower and conrows.upper arguments to specify the bounds on 
	the constraints. \code{conrows} is a matrix of npar(object) columns and
	one row for each constraint (ie a vector in the case of a single
	constraint). 
	
	\code{llratio} performs a log-likelihood ratio test on two
	\code{fit}ted models; the first object should have the largest degrees 
	of freedom. 
	
}

\value{
	
	\code{fit} returns an object of class \code{\link[depmix.fitted-class]{depmix.fitted}} which contains 
	the original depmix object, and further has slots: 
	
	\item{message}{: Convergence information.}
	
	\item{conMat}{: The constraint matrix A, see details.}

	\item{posterior}{: The posterior state probabilities.}
	
	The print method shows the \code{message} and the summary method
	shows the parameter estimates.
	
}

\examples{

data(speed)
# 2-state model on the RTs of the speed data with random 
# starting values for the transition pars (without those EM does not get off the ground)
set.seed(1)
mod <- depmix(rt~1,data=speed,nstates=2,trstart=runif(4))
# fit the model
mod1 <- fit(mod)
mod1 # to see the logLik and optimization information
# to see the parameters
summary(mod1)

}

\author{Ingmar Visser}

\keyword{methods}
