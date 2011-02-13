\name{summary.gibbs}
\alias{summary.gibbs}
\alias{print.summary.gibbs}


\title{Summary of results from the Gibbs Sampler}
\description{
  This function is a \code{method} for class \code{gibbs}. It provides a summary of the results from the Gibbs sampler - descriptive statistics of
parameters and functions of parameters, and posterior predictive checks.
}

\usage{
\method{summary}{gibbs}(object, start, end, D, CORR=TRUE, PROBT1=TRUE, 
DIC=TRUE, PROBD=TRUE,  probs= c(0.025, 0.50, 0.975), \dots)
}


\arguments{
  \item{object}{an object of class \code{gibbs}, a result of a call to \code{gs}.}  
  \item{start}{the Gibbs sampler iteration from which to start using
  samples for the summary.}  
 \item{end}{the Gibbs sampler iteration from which to stop using
  samples for the summary.}  
  \item{D}{a list defining latent variables with respect to which we want to calculate sensitivity and specificity. A latent variable
could be defined by a single latent class or by two or more latent classes. See example for details. The default value is
a latent variable that is defined by the first latent class.}
  \item{probs}{numeric vector of probabilities with values in [0,1], indicating quantiles to be reported as summary statistics. 
  The default vector is c(0.025, 0.50, 0.975).}   
  \item{CORR}{logical. Set to TRUE to calculate the observed and expected correlation, and the correlation residuals. Otherwise set to FALSE.}
  \item{PROBT1}{logical. Set to TRUE to calculate P(positive result|latent class) for each variable (test). Otherwise set to FALSE.}
  \item{DIC}{logical. Set to TRUE to calculate the log-likelihood, the DIC, pD and effective number of parameters. Otherwise set to FALSE.}
  \item{PROBD}{logical. Set to TRUE to calculate P(latent class|combination of observed variables. Otherwise set to FALSE.}    
  \item{\dots}{further arguments passed to or from other methods.}  
}

 

\value{  
  \item{iter}{the total number of iterations of the Gibbs sampler on which the results are based.}
  \item{corr}{a list of 'corr' and 'pCorr', giving the observed and predicted correlation of all the pairs of tests for the iterations from 'burnin+1' to 'iter'.}
   \item{probT1}{an array giving P(positive result|latent class) for each observed variable, based on iterations from 'burnin+1' to 'iter'.}
  \item{probD}{an array giving P(latent class|combination of observed variables), based on iterations from 'burnin+1' to 'iter'.}
  \item{ss}{a list giving the sensitivities and specificities with respect to the variables specified in argument \eqn{D}.}  
  \item{Dbar}{posterior mean deviance.}
  \item{pD}{effective number of parameters.}
  \item{DIC}{Deviance Information Criterion (DIC) = Dbar + pD.}  
    \item{summary.prev}{a matrix of nClass rows and (length(probs)+1) columns. This matrix give the quantiles corresponding to the probabilities in the vector 
'probs' and the means for the probabilities for each latent class.}
  \item{summary.a}{a matrix of 'nA' rows and (length(probs)+1) columns. This matrix give the quantiles and the means for the \eqn{a} parameters.}
  \item{summary.b}{a matrix of 'nB' rows and (length(probs)+1) columns. This matrix give the quantiles and the means for the \eqn{b} parameters.}
  \item{summary.prEgO}{a matrix of 'nClass' rows with entries Pr(E>O), the probability that the predicted percentage agreement is greater than the observed agreement for every pair of variables, within each latent class.}
  \item{summary.agree}{a list of 'agree' and 'pAgree', the observed and predicted percentage agreements for every pair of variables.}
  \item{summary.corr}{a list of  'corr', 'pCorr' and 'resid', the observed and predicted correlation and their difference for every pair of variables.}
  \item{summary.pred}{a matrix giving the predicted frequency of each combination of variables.}
  \item{summary.probD}{a list of length nClass giving P(latent class|combination of observed variables) for all observed combinations of variables. }
  \item{summary.corr}{a list of 'corr','pCorr' and 'resid',
  the observed and predicted correlation and their difference for every 
  pair of variables.}
 \item{summary.probT1}{a list of 'Class1', ..., 'ClassN' giving Pr(positive result|latent class), where N is the number of the latent classes. }
}
\seealso{ \code{\link{gs}}}
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
summary(tlcm.2re,5)
}
\keyword{methods}

