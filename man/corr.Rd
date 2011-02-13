\name{corr}
\alias{corr}

\title{calculating the correlation residuals}
\description{
  This function  calculates the observed and predicted correlation residuals.
}
\usage{
corr(x, start, end,\dots)
}


\arguments{
 \item{x}{an object of class \code{gibbs}, a result of a call to \code{gs}.}  
 \item{start}{the number of iterations from which to start using
 samples to calculate  the correlation residuals.}  
  \item{end}{the number of iterations from which to stop using
 samples to calculate  the correlation residuals.}  
 \item{\dots}{further arguments passed to or from other methods.}  
}


\value{
An list of \code{corr} and \code{pCorr} for observed and predicted correlation residuals, 
respectively. 
}


\seealso{ \code{\link{gs}}  }
\examples{
##--------------- hlcm basic ---------------------------
cstA <- matrix(0, 3,4)                 
cstA[1,] <- c(1,2,3,4)                                       
cstA[2,] <- c(1,2,5,6)                                     
cstA[3,] <- c(7,8,5,6)    
cstB     <- matrix(0, 3,4)                                
init2    <- list(prev=c(0.2, 0.5, 0.3), a= c(rep(2, 4), rep(-2, 4)))
lcmr.obj <- lcmr(sample.Chlamydia, cstA, cstB)
hlcm.basic <- gs(lcmr.obj, iter=10,  init=init2)
re.corr    <- corr(hlcm.basic, 6)
}
\keyword{methods}

