\name{mu}
\alias{mu}

\title{calculating the probability of a positive result of a test}
\description{
  This function  calculates the probability of a positive result of a test conditional
 on the latent class  and random effects.
}
\usage{
mu(x, start, end, \dots)
}


\arguments{
 \item{x}{an object of class \code{gibbs}, a result of a call to \code{gs}.}  
 \item{start}{the Gibbs sampler iteration from which to 
  start using samples for  the calculation.}  
   \item{end}{the Gibbs sampler iteration from which to 
 stop using samples for  the calculation.}  
 \item{\dots}{further arguments passed to or from other methods.}  
}


\value{
A list.  Its  length is equal to the  number of interations after the burn-in period. 
Each element of this list is array of \code{N*nTest*nClass}.
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
re.mu <- mu(hlcm.basic, 6)
}
\keyword{methods}

