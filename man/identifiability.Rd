\name{identifiability}
\alias{identifiability}

\title{Check the identifiability}

\description{  
The function \code{identifiability}  is used to check the identifiability of a \code{lcmr} model by investigating
whether the necessary and sufficient conditions for identifiability are satisfied.  The necessary condition is:
number of parameters <= number of degrees of freedom. The sufficient condition is that the Jacobian matrix has full rank.   
}

\usage{
identifiability(lcmrObj)
}

\arguments{  
  \item{lcmrObj}{an object of class \code{lcmr}.}
}


\value{
   An list is returned containing  the following components:
   \item{nPar:}{the number of parameters.}   
   \item{ndf:}{the number of degrees of freedom.}
   \item{rank:}{the rank of the Jacobian matrix.}
   \item{identifiable:}{a logical value. 0: the necessary condition is not satisfied; 1: both the necessary and sufficient conditions are satisfied;
    2: the necessary condition is satisfied, but the sufficient condition is not satisfied.}   
}

\references{
Anton K.Formann (2003), Latent Class Model Diagnosis from a Frequentist Point of View, \code{Biometrics}, 59, 189-196

Leo A. Goodman (1974), Exploratory Latent Structure Analysis Using Both Identifiable and Unidentifiable Models, \code{Biometrika}, 61, 215-231
}

\examples{
init <- list(prev=c(0.2, 0.8), a= c(rep(2, 4), rep(-2, 4)), b=c(1, 1))
cstB <- matrix(0, 2,4)
cstB[1, ] <-  c(1, 1, 2, 2)    
cstB[2, ] <-  c(0, 0, 0, 0)    

cstA <- matrix(0, 2,4)
cstA[1,] <- c(1 ,2 ,3 ,4)
cstA[2,] <- c(5, 6, 7, 8)

lcmr.obj <- lcmr(sample.Chlamydia, cstA, cstB)
identifiability(lcmr.obj)    
}
\keyword{methods}

