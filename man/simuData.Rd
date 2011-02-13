\name{simuData}
\alias{simuData}

\title{Generate simulated data set}
\description{
  This function can be used to simulate data from a latent class model assuming conditional independence between the observed dichotomous variables.
}
 

\usage{
simuData(N, prev, trueP) 
}

\arguments{    
  \item{N}{sample size.}  
  \item{prev}{a vector of prevalence for each class.}  
  \item{trueP}{a matrix of nClass by nTest giving the probability of a positive result on each variable within each latent class P(Variable=1|Latent class), where nTest is the number of observed dichotomous variables.}
}
 

\value{  
\item{dat}{the simulated data set, a matrix with (nTest+1) columns, one column for each variable and the last column for the frequency of each combination of results.}
\item{prev}{a vector giving the proportion in each latent class in the simulated dataset.}
\item{p}{a vector giving estimates of the probabilities in \eqn{trueP} based on the simulated dataset.}
}


\examples{
N <- 1000;  # the number of subjects
nClass <- 2;     # the number of classes            
prev   <- c(0.95, 0.05)    # prevalence of 1: D+; 2: D-     
P <- 6;                  # the number of tests
trueP <- matrix(c(.85,.85,.85,.75,.75,1,.05,.05,.05,.01,.01,0), nClass, P, byrow=TRUE)
simu.data <- simuData(N,  prev,trueP)
simu.data$dat             # the simulated dataset}
\keyword{methods}

