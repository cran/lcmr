\name{gs}
\alias{gs}
\alias{update.gibbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Gibbs sampler}
\description{
  \code{gs} is used to estimate the latent class model 
  using a Bayesian approach via a Gibbs sampler. 
  This function \code{update} is a \code{method} for class \code{gibbs} 
  objects, which is used  to update the  object that we get from the 
  function \code{gs} by running more iterations.
  } 


\details{
  The probability of a positive result on each observed variable takes the form pnorm(\eqn{a+b*r}), where \eqn{r} is a N(0,1) random effect and \eqn{a}
  and \eqn{b} are unknown parameters. The \eqn{a} parameters follow \eqn{Normal(prior.mA, prior.sdA)} prior distributions. The 
  \eqn{b} parameters follow \eqn{Uniform(priorB1, priorB2)} prior distributions. The multinomial probability of the different
latent classes follows a \eqn{Dirichlet(priorPrev[1], ..., priorPrev[nClass])} prior distribution. 
}

\usage{
gs(lcmr.obj, iter,  init, prior.sdA, prior.mA, priorPrev, 
 priorA=list(l=-1000, u=1000), priorB=list(l=0, u=5), \dots) 
\method{update}{gibbs}(object, newiter, \dots) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{  
  \item{lcmr.obj}{an object of class \code{lcmr}.}
  \item{iter}{the number of iterations of the Gibbs sampler.}
  \item{init}{a list of initial values for the parameters in the model, i.e. the probability of each latent class, a's and b's.}
  \item{prior.sdA}{a vector with the prior standard deviations of the \eqn{a}'s. The default value is a vector of 1's. The length
of the vector is equal to the number of unique \eqn{a} parameters specified in the lcmr.obj.}
  \item{prior.mA}{a vector with the prior means of the \eqn{a}'s. The default value is a vector of 0's. The length
of the vector is equal to the number of unique \eqn{a} parameters specified in the lcmr.obj.}
  \item{priorPrev}{a vector with the parameters of the Dirichlet prior distribution for the probability of each latent classes. The default value will be a vector of all 1's. The length of the vector is equal to the number of latent classes.}
  \item{priorA}{a list with the lower limit and the upper limit of the range of the \eqn{a}'s. 
\eqn{priorA$l} and \eqn{priorA$u} are either  two vectors of length equal to the number of \eqn{a}'s or two scalars 
when  all a's have the same range. The default value is [-1000, 1000].}
  \item{priorB}{a list with the lower limit and the upper limit of the range of the \eqn{b}'s. 
\eqn{priorB$l} and \eqn{priorB$u} are either  two vectors of length equal to the number of \eqn{a}'s or two scalars 
when  all a's have the same range. The default value is [0, 5].}
 \item{object}{an object of class \code{gibbs}, a result of a call to \code{gs}.}  
  \item{newiter}{the number of new iterations.}  
\item{\dots}{further arguments passed to or from other methods.}
}

\value{
An object of class \code{gibbs}. A list containing at least the following components:
   \item{N}{the number of patients.}   
   \item{nTest}{the number of dichotomus variables.}
   \item{nCell}{the number of observed combinations of the dichotomous variables.}
   \item{nClass}{the number of the latent classes.}
   \item{nA}{the number of \eqn{a} parameters. }
   \item{nB}{the number of \eqn{b} parameters. }
   \item{nR}{the number of random effects. }
   \item{allCell}{the frequency of each possible combination of the dichotomous variables, ranging from all 0's to all 1's.}
  \item{data}{the \eqn{data} argument}
  \item{DAT}{a matrix of nTest columns and N rows. The argument \eqn{data} transformed into a format where each line corresponds to a single observation.}
  \item{cstA}{the \eqn{cstA} argument.}
  \item{cstB}{the \eqn{cstB} argument}
\item{iter}{number of iteration of the Gibbs sampler.}
\item{burnin}{number of burnin iterations. }
\item{init}{initial values used to commence the Gibbs sampler}
\item{prev}{a matrix of 'iter' rows and 'nClass' columns giving  the probability for each latent class for all the iterations.}
\item{a}{a matrix of 'iter' rows and 'nA' columns giving  the sample of \eqn{a}'s drawn by the Gibbs sampler.}
\item{b}{a matrix of 'iter' rows and 'nB' columns giving  the sample of \eqn{b}'s drawn by the Gibbs sampler.}
\item{agree}{a list of 'pAgree' and 'agree', giving the observed and predicted percentage agreements of all the pairs of tests.}
\item{predmat}{a matrix giving the predicted frequency of each combination of the observed variables.}
\item{logLik}{log-likelihood for all the iterations.}
}
  
\examples{
##---------- tlcm basic --------------------------
cstA <- matrix(0, 2, 4)
cstA[1,] <- c(1 ,2 ,3 ,4)
cstA[2,] <- c(5, 6, 7, 8)
cstB <- matrix(0, 2, 4)
lcmr.obj <- lcmr(sample.Chlamydia, cstA, cstB)
init <- list(prev=c(0.2, 0.8), a= c(rep(2, 4), rep(-2, 4)))
tlcm.basic  <- gs(lcmr.obj,  iter=10,  init=init)                          

##---------- tlcm.2re--------------------------
cstA <- matrix(0, 2,4)
cstA[1,] <- c(1 ,2 ,3 ,4)
cstA[2,] <- c(5, 6, 7, 8)

cstB1 <- matrix(0, 2, 4)
cstB1[1, ] <-  c(1, 1, 0, 0)    
cstB2 <- matrix(0, 2, 4)
cstB2[1, ] <-  c(0, 0, 2, 2)   
cstB <- list(B1=cstB1, B2=cstB2)
init <- list(prev=c(0.2, 0.8), a= c(rep(2, 4), rep(-2, 4)), b=c(1, 1))
lcmr.obj <- lcmr(sample.Chlamydia, cstA, cstB)
tlcm.2re  <- gs(lcmr.obj,  iter=10,  init=init)               

##--------------- hlcm basic ---------------------------
cstA <- matrix(0, 3,4)                 
cstA[1,] <- c(1,2,3,4)                                       
cstA[2,] <- c(1,2,5,6)                                     
cstA[3,] <- c(7,8,5,6)    
cstB <- matrix(0, 3,4)                                
init2 <- list(prev=c(0.2, 0.5, 0.3), a= c(rep(2, 4), rep(-2, 4)))
lcmr.obj <- lcmr(sample.Chlamydia, cstA, cstB)
hlcm.basic <- gs(lcmr.obj, iter=10, init=init2)

#---------------- hlcm+re 3b --------------------------
cstA <- matrix(0, 3,4)                 
cstA[1,] <- c(1,2,3,4)                                       
cstA[2,] <- c(1,2,5,6)                                     
cstA[3,] <- c(7,8,5,6)   
cstB <- matrix(0, 3, 4)                 
cstB[1, ] <- c(1, 1, 0, 0)             
cstB[2, ] <- c(2, 2, 0, 0)             
cstB[3, ] <- c(3, 3, 0, 0)                                 
init2 <- list(prev=c(0.2, 0.5, 0.3), a= c(rep(2, 4), rep(-2, 4)), b=c(1, 1, 1))
lcmr.obj   <- lcmr(sample.Chlamydia, cstA, cstB)
hlcm.re.3b <- gs(lcmr.obj, iter=10,  init=init2)

##--------------- hlcm+re 2b ----------------------------
cstA <- matrix(0, 3,4)                 
cstA[1,] <- c(1,2,3,4)                                       
cstA[2,] <- c(1,2,5,6)                                     
cstA[3,] <- c(7,8,5,6) 

cstB <- matrix(0, 3,4)                 
cstB[1, ] <- c(1, 1, 0, 0)             
cstB[2, ] <- c(2, 2, 0, 0)             
cstB[3, ] <- c(0, 0, 0, 0)             
                                                              
init2 <- list(prev=c(0.2, 0.5, 0.3), a= c(rep(2, 4), rep(-2, 4)), b=c(1, 1))
lcmr.obj   <- lcmr(sample.Chlamydia, cstA, cstB)
hlcm.re.2b <- gs(lcmr.obj, iter=10,  init=init2)

##---------------------- hlcm+re 1b -----------------------
cstA <- matrix(0, 3,4)                 
cstA[1,] <- c(1,2,3,4)                                       
cstA[2,] <- c(1,2,5,6)                                     
cstA[3,] <- c(7,8,5,6)   

cstB <- matrix(0, 3,4)                 
cstB[1, ] <- c(1, 1, 0, 0)             
cstB[2, ] <- c(0, 0, 0, 0)             
cstB[3, ] <- c(0, 0, 0, 0)                             
init2 <- list(prev=c(0.2, 0.5, 0.3), a= c(rep(2, 4), rep(-2, 4)), b=c(1))
lcmr.obj   <- lcmr(sample.Chlamydia, cstA, cstB)
hlcm.re.1b <- gs(lcmr.obj, iter=10,  init=init2)
##--------------------- update hlcm.basic -----------------------
update.hlcm.basic <- update(hlcm.basic, 10)
}

\keyword{methods}

