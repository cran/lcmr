\name{logLik.gibbs}
\alias{logLik.gibbs}

\title{The log likelihood function and the deviance information criterion (DIC)}
\description{
  This function is to return the values of the log likelihood function and the deviance information criterion (DIC).
}
\usage{
 \method{logLik}{gibbs}(object, start, end, \dots)
}


\arguments{
 \item{object}{an object of class \code{gibbs}, a result of a call to \code{gs}.}  
 \item{start}{the Gibbs sampler iteration from which to 
  start returning the log likelihood values.}  
   \item{end}{the Gibbs sampler iteration from which to 
  stop returning the log likelihood values.}  
 \item{\dots}{further arguments passed to or from other methods.}  
}


\value{
A list of 
\item{logLike}{A vector of log-likelihood values from the iteration number  \code{burnin+1}. }
\item{sumy.logLik}{a summary of the values of log-likelihood function.}
\item{dbar}{the average of deviance \code{D(theta)} over the samples of parameter  \code{theta}.}
\item{dhat}{the value of deviance \code{D} evaluated at the average of the samples of \code{theta}.}
\item{pd}{the effective number of parameters.}
\item{dic}{the value of deviance information criterion.}
}


\seealso{ \code{\link{gs}}  }
\examples{
\dontshow{
##--------------- hlcm basic ---------------------------
cstA <- matrix(0, 3,4)                 
cstA[1,] <- c(1,2,3,4)                                       
cstA[2,] <- c(1,2,5,6)                                     
cstA[3,] <- c(7,8,5,6)    
cstB <- matrix(0, 3,4)                                
init2 <- list(prev=c(0.2, 0.5, 0.3), a= c(rep(2, 4), rep(-2, 4)))
lcmr.obj <- lcmr(sample.Chlamydia, cstA, cstB)
hlcm.basic <- gs(lcmr.obj, iter=10, init=init2)
}
## --- Continuing the Example from  '? gs':
re.logLik <- logLik(hlcm.basic, 6)
}
\keyword{methods}

