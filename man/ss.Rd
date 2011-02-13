\name{ss}
\alias{ss}

\title{Calculating sensitivity and specificity}
\description{
  This function calculates the sensitivity and specificity for
  the  object obtained from the function \code{gs}.
}
 

\usage{
ss(x, start, end, D=1, probs= c(0.025, 0.50, 0.975), \dots)
}


\arguments{
  \item{x}{an object of class \code{gibbs}, a result of a call to \code{gs}.}  
  \item{start}{the Gibbs sampler iteration from which to 
  start using samples for calculating the sensitivity and specificity.}  
  \item{end}{the Gibbs sampler iteration from which to 
  stop using samples for calculating the sensitivity and specificity.}  
  \item{D}{a vector of latent class subsets for which we want to calculate the sensitivities and specificities. For example,
   in the 3-classes model in the following reference, we can indicate D=1 or D= c(1, 2).}
  \item{probs}{numeric vector of probabilities with values in [0,1], indicating quantiles to be reported as summary statistics. 
  The default vector is c(0.025, 0.50, 0.975).}     
  \item{\dots}{further arguments passed to or from other methods.}  
}

 

\value{  
\item{sensitivity}{Sensitivities.}
\item{specificity}{Specificities.}  
}


\examples{
\dontshow{
######## tlcm.2re
init <- list(prev=c(0.2, 0.8), a= c(rep(2, 4), rep(-2, 4)), b=c(1, 1))
cstB <- matrix(0, 2,4)
cstB[1, ] <-  c(1, 1, 2, 2)    
cstB[2, ] <-  c(0, 0, 0, 0)    

cstA <- matrix(0, 2,4)
cstA[1,] <- c(1 ,2 ,3 ,4)
cstA[2,] <- c(5, 6, 7, 8)
cstR <- c(1, 1, 2, 2)
lcmr.obj <- lcmr(sample.Chlamydia, cstA, cstB)
tlcm.2re <- gs(lcmr.obj, iter=10, init=init)
}
## --- Continuing the Example from  '? gs':
re.ss <- ss(tlcm.2re, 6)
}
\keyword{methods}

