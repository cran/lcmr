\name{lcmr}
\alias{lcmr}

\title{Create a lcmr object}
\description{
  \code{lcmr} is a function to create an object of class \code{lcmr}, which represents a latent class model with random effects. It 
  is used as an argument in the function \code{gs}.}

\details{
  The probability of a positive result on each observed variable in the latent class model takes the form pnorm(\eqn{a+b*r}), where \eqn{r} is a N(0,1) random effect and \eqn{a}
  and \eqn{b} are unknown parameters.
}

\usage{
lcmr(data, cstA, cstB) 
}

\arguments{  
  \item{data}{a matrix with (nTest+1) columns, one column for each variable and the last column for the frequency of each combination of results.}
  \item{cstA}{a matrix of nClass rows and nTest columns specifying the number of unique the \eqn{a} parameters. nClass is the number of latent classes in the model
 and nTest is the number of observed dichotomous varables. Each unique parameter is denoted by a different positive integer in this matrix.}
  \item{cstB}{a list in which every element is a matrix of nClass rows and nTest columns. The number of matrices in the list is equal to the number of random effects
in the latent class model. The elements of the matrices are used to specify the number of unique \eqn{b} parameters. Each unique parameter is denoted by a positive integer. 
If the [i,j]th entry in a matrix is 0 it implies that the jth variable is not influenced by the random effect in the ith latent class. If the number of random effects is 0, 
cstB is a matrix of 0's of nClass rows and nTest columns.} 
}


\value{
   An object of class \code{lcmr} is a list containing at least the following components:
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
}


\examples{

## 4 tests, two latent classes, 2 random effects 
cstB <- matrix(0, 2,4)
cstB[1, ] <-  c(1, 1, 2, 2)    
cstB[2, ] <-  c(0, 0, 0, 0)    
cstA <- matrix(0, 2,4)
cstA[1,] <- c(1 ,2 ,3 ,4)
cstA[2,] <- c(5, 6, 7, 8)
lcmr.obj <- lcmr(sample.Chlamydia, cstA, cstB)

## 4 tests, 3 latent classes, no contraints on a's, no random effects 
cstA <- matrix(0, 3,4)                 
cstA[1,] <- c(1,2,3,4)                                       
cstA[2,] <- c(5,6,7,8)                                     
cstA[3,] <- c(9,10,11,12)    
cstB <- matrix(0, 3,4)                                
lcmr.obj <- lcmr(sample.Chlamydia, cstA, cstB)

##--------------- 4 tests, 3 latent classes, contraints on a's, no random effects ---------------------------
cstA <- matrix(0, 3,4)                 
cstA[1,] <- c(1,2,3,4)                                       
cstA[2,] <- c(1,2,5,6)                                     
cstA[3,] <- c(7,8,5,6)    
cstB <- matrix(0, 3,4)                                
lcmr.obj <- lcmr(sample.Chlamydia, cstA, cstB)

}
\keyword{methods}

