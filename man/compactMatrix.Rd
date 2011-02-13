\name{compactMatrix}
\alias{compactMatrix}

\title{ Make the data set compact, or extract a sub data set}
\description{
  This function can be used to extract the cross-tabulation between a subset of the dichotomous variables in the larger data set.
} 

\usage{
compactMatrix(M, cols) 
}

\arguments{
  \item{M}{a matrix with (nTest+1) columns, one column for each variable and the last column for the frequency of each combination of results.}  
  \item{cols}{a vector indicating the columns numbers, ranging from 1 to nTest, that are to be extracted .}    
  
}
 

\value{  
The subset of the data set indicated by the selected columns. 
}

\examples{
\dontshow{
N <- 1000;  # the number of subjects
nClass <- 2;     # the number of classes            
prev <- c(0.95, 0.05)    # prevalence of 1: D+; 2: D-     
P <- 6;                  # the number of tests
trueP <- matrix(c(.85,.85,.85,.75,.75,1,.05,.05,.05,.01,.01,0), nClass, P, byrow=TRUE)
simu.data <- simuData(N,  prev,trueP) 
}
## --- Continuing the Example from  '? simuData':
# Create data sets with (T1, T2, T3, T4), (T1, T2, T4, T5) 
# and (T1, T4, T5, T6) from the data set, simu.data$dat
sub1 <- compactMatrix(simu.data$dat, c(1, 2, 3, 4)) 		# (T1, T2, T3, T4)
sub2 <- compactMatrix(simu.data$dat, c(1, 2, 4, 5))       	# (T1, T2, T4, T5)
sub3 <- compactMatrix(simu.data$dat, c(1, 4, 5, 6))			# (T1, T4, T5, T6)
}
\keyword{methods}

