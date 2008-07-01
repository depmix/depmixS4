\name{simulate}

\docType{method}

\alias{simulate}

\alias{simulate,response-method}
\alias{simulate,depmix-method}
\alias{simulate,GLMresponse-method}
\alias{simulate,transInit-method}

\alias{simulate,BINOMresponse-method}
\alias{simulate,GAMMAresponse-method}
\alias{simulate,MULTINOMresponse-method}
\alias{simulate,MVNresponse-method}
\alias{simulate,NORMresponse-method}
\alias{simulate,POISSONresponse-method}

\title{Methods to simulate from depmix models}

\description{

	Random draws from \code{depmix} objects.

}

\usage{

  \S4method{simulate}{depmix}(object, nsim=1, seed=NULL, ...)
    
	\S4method{simulate}{response}(object, nsim=1, seed=NULL, times, ...)
	
	\S4method{simulate}{GLMresponse}(object, nsim=1, seed=NULL, times, ...)
	
	\S4method{simulate}{transInit}(object, nsim=1, seed=NULL, times, is.prior=FALSE, ...)
}

\arguments{
	
	\item{object}{Object to generate random draws.  An object of class
	\code{depmix}, \code{response}, \code{transInit}}

	\item{nsim}{The number of draws (one draw simulates a data set of the
	size that is defined by ntimes); defaults to 1.}
	
	\item{seed}{Set the seed.}
	
	\item{times}{(optional) An indicator vector indicating for which times in the
	complete series to generate the data.  For internal use.}
	
	\item{is.prior}{For \code{transInit} objects, indicates whether it is a prior 
  (init) model, or not (i.e., it is a transition model)}
	
	\item{...}{Not used currently.}
}

\details{

	For a \code{depmix} model, simulate generates \code{nsim} random state
	sequences, each of the same length as the observation sequence in the
	\code{depmix} model (i.e., \code{sum(ntimes(object))}.  The state
	sequences are then used to generate \code{nsim} observation sequences
	of thee same length.

	Setting the \code{times} option selects the time points in the total
	state/observation sequence (i.e., counting is continued over ntimes).
	Direct calls of simulate with the \code{times} option are not recommended.

}

\value{

	For a \code{depmix} object, a new object of class \code{depmix.sim}. 
	
	For a \code{transInit} object, a state sequence.
	
	For a \code{response} object, an observation sequence.

}

\author{Maarten Speekenbrink}

\keyword{methods}