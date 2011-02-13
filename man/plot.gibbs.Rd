\name{plot.gibbs}
\alias{plot.gibbs}

\title{plotting a Gibbs sampler}
\description{
  This function is a \code{method} for class \code{gibbs}  objects.
  It can be used to obtain trace plots that serve as a crude check of convergence. It
can also be used to obtain plots corresponding to the posterior predictive check for
conditional dependence between pairs of tests.
}
\usage{
\method{plot}{gibbs}(x, option=1, ask=TRUE, \dots)
}


\arguments{
 \item{x}{an object of class \code{gibbs}, a result of a call to \code{gs}.}  
 \item{option}{option=1: trace plots; option=2: posterior predictive check.}  
 \item{ask}{TRUE: ask before displaying the next plot; FALSE: do not ask.}
 \item{\dots}{further arguments passed to or from other methods.}  
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
plot(hlcm.basic)
}
\keyword{methods}

