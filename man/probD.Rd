\name{probD}
\alias{probD}

\title{calculating the Pr(Latent class|Pattern)}
\description{
  This function  returns the probability of a latent class given a pattern: Pr(Latent class|Pattern).
}
\usage{
probD(x, start, end, probs=c(0.025, 0.50, 0.975), \dots)
}


\arguments{
 \item{x}{an object of class \code{gibbs}, a result of a call to \code{gs}.}  
 \item{start}{the Gibbs sampler iteration from which to 
  start using samples for calculating 
 Pr(Latent class|Pattern).}  
 \item{end}{the Gibbs sampler iteration from which to 
  stop using samples for calculating 
 Pr(Latent class|Pattern).} 
 \item{probs}{numeric vector of probabilities with values in [0,1], indicating quantiles to be reported as summary statistics. 
  The default vector is c(0.025, 0.50, 0.975).}   
 \item{\dots}{further arguments passed to or from other methods.}  
}


\value{
An array of (number of patterns)*nClass*(number of iterations). 
}


\seealso{ \code{\link{gs}}  }
\examples{
\dontshow{
##--------------- hlcm basic ---------------------------
cstA <- matrix(0, 3,4)                 
cstA[1,] <- c(1,2,3,4)                                       
cstA[2,] <- c(1,2,5,6)                                     
cstA[3,] <- c(7,8,5,6)    
cstB     <- matrix(0, 3,4)                                
init2    <- list(prev=c(0.2, 0.5, 0.3), a= c(rep(2, 4), rep(-2, 4)))
lcmr.obj <- lcmr(sample.Chlamydia, cstA, cstB)
hlcm.basic <- gs(lcmr.obj, iter=10, init=init2)
}
## --- Continuing the Example from  '? gs':
re.probD <- probD(hlcm.basic, 6)
}
\keyword{methods}

