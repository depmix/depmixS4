
\name{depmix}

\docType{methods}

\alias{depmix}
\alias{depmix,ANY-method}

\alias{show}
\alias{summary}

\alias{show,depmix-method}
\alias{summary,depmix-method}

\title{Dependent Mixture Model Specifiction}

\description{
	
	\code{depmix} creates an object of class \code{depmix}, a dependent
	mixture model, otherwise known as hidden Markov model.  For a short
	description of the package see \code{\link{depmixS4}}.
	
}

\usage{
	
	depmix(response, data=NULL, nstates, transition=~1, family=gaussian(), 
		prior=~1, initdata=NULL, respstart=NULL, trstart=NULL, instart=NULL,
		ntimes=NULL,...)	
	
}

\arguments{
	
	\item{response}{The response to be modeled; either a formula or a list 
		of formulae in the multivariate case; this interfaces to the glm
		distributions. See 'Details'.}
	
	\item{data}{An optional data.frame to interpret the variables in
		the response and transition arguments.}
	
	\item{nstates}{The number of states of the model.}

	\item{transition}{A one-sided formula specifying the model for the 
		transitions. See 'Details'.}
	
	\item{family}{A family argument for the response. This must be a list
		of family's if the response is multivariate.}
	
	\item{prior}{A one-sided formula specifying the density for the prior
		or initial state probabilities.}
	
	\item{initdata}{An optional data.frame to interpret the variables
		occuring in prior. The number of rows of this data.frame must be
		equal to the number of cases being modeled. See 'Details'.}
	
	\item{respstart}{Starting values for the parameters of the response
		models.}
	
	\item{trstart}{Starting values for the parameters of the transition
		models.}
	
	\item{instart}{Starting values for the parameters of the prior or
		initial state probability model.}

	\item{ntimes}{A vector specifying the lengths of individual, ie
		independent, time series. If not specified, the responses are
		assumed to form a single time series. If the data argument has an
		attribute ntimes, then this is used.}
		
	\item{...}{Not used currently.}
		
}

\details{ 

	The function \code{depmix} creates an S4 object of class \code{depmix},
	which needs to be fitted using \code{\link{fit}} to optimize the
	parameters.
	
	The response model(s) are created by call(s) to \code{\link{response}}
	providing the response formula and the family specifying the error
	distribution.  If response is a list of formulae, the \code{response}'s
	are assumed to be independent conditional on the latent state.
	
	The transitions are modeled as a multinomial logistic model for each
	state.  Hence, the transition matrix can be modeled as time-dependent,
	depending on predictors.  The prior density is also modeled as a
	multinomial logistic.  Both are created by calls to
	\code{\link{transInit}}.
	
	Starting values may be provided by the respective arguments.  The order
	in which parameters must be provided can be easily studied by using the
	\code{\link{setpars}} function.
	
	Linear constraints on parameters can be provided as argument to the
	\code{\link{fit}} function.
	
	The print function prints the formulae for the response, transition and
	prior models along with their parameter values.
	
}

\value{
	
\code{depmix} returns an object of class \code{depmix} which has the
following slots:

	\item{response}{A list of a list of response models; the first
		index runs over states; the second index runs over the independent 
		responses in case a multivariate response is provided.}
	
	\item{transition}{A list of \code{transInit} models, ie multinomial
		logistic models with length the number of states.}

	\item{prior}{A multinomial logistic model for the initial state
		probabilities.}

	\item{dens,trDens,init}{See \code{\link{depmix-class}} help for
		details.  For internal use.}
	
	\item{stationary}{Logical indicating whether the transitions are
		time-dependent or not; for internal use.}
	
	\item{ntimes}{A vector containing the lengths of independent time
		series; if data is provided, sum(ntimes) must be equal to
		nrow(data).}
	
	\item{nstates}{The number of states of the model.}
	
	\item{nresp}{The number of independent responses.}

	\item{npars}{The total number of parameters of the model.  Note: this
	is \emph{not} the degrees of freedom because there are redundancies in the
	parameters, in particular in the multinomial models for the transitions
	and prior probabilities.}

}

\author{Ingmar Visser & Maarten Speekenbrink}

\seealso{

	\code{\link{fit}}, \code{\link{transInit}}, \code{\link{response}},
	\code{\link{depmix-methods}} for accessor functions to \code{depmix}
	objects.

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
data(speed)
mod <- depmix(list(rt~1,corr~1),data=speed,nstates=2,family=list(gaussian(),multinomial()))
# print the model, formulae and parameter values
mod

}

\keyword{methods}


