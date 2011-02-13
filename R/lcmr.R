###            Fit a latent class model with random effects
###
### Copyright 2005-2011  Liangliang Wang <l.wang@stat.ubc.ca>,
###                      Nandini Dendukuri <nandini.dendukuri@mcgill.ca>.
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#========================lcmr===================================
# data: 1, 2, ... nTest; count

lcmr   <- 
## latent class model with random effects 
    function(data, cstA, cstB)
{		       
  #----------- check the input arguments---------------
  #----------- check 'data'
  if(is.null(dim(data))) stop("Please check your input data set\n")                           
  nTest <- dim(data)[2]-1   # number of tests 	
  data <- compactMatrix(data, 1:nTest)
  if(nTest < 2) stop("The number of tests should be greater than 1!\n")        
  nCell <- dim(data)[1]     # number of patterns or cells
  if(nCell < 2) stop("The number of patterns should be greater than 1!\n") 
  if(sum(data[, 1:nTest]!=0 &data[, 1:nTest]!=1)>0)
     stop("Results of tests should be either 0 or 1!\n") 
  if(sum(data[, nTest+1]<0 )>0)
     stop("The number of subjects for each pattern should be positive!\n")
  #--------- check 'cstA'
  if(class(cstA)!="matrix") stop("cstA is a matrix of nClass by nTest!\n")
  nA <- length(unique(cstA[cstA > 0])) 
  if(nA==0) stop("The number of parameter a's should not be 0!\n")
  #--------- check 'cstB'
  if(missing(cstB)){
    stop("cstB is a list in which every element is a matrix of nClass*nTest. The length of cstB   is    equal to the number of random effects. If the number of random effects is 1, cstB can bea matrix of nClass*nTest.\n")	
    }  
  if(class(cstB)!="list") {
    if(class(cstB)=="matrix" & all(dim(cstB)==dim(cstA))){ cstB = list(B1=cstB)}
    else{stop("cstB is a list in which every element is a matrix of nClass*nTest. The length of cstB   is    equal to the number of random effects. If the number of random effects is 1, cstB can bea matrix of nClass*nTest.\n")} 
  }
  #----------- end of checking the input arguments---------------
	
  test.names <- colnames(data)[1:nTest]
  N          <- sum(data[, nTest+1])    # total number of observations	
  nClass     <- nrow(cstA)
  nR         <- length(cstB)
  unlistcstB <- unlist(cstB)
  nB         <- length(unique(unlistcstB[unlistcstB > 0]));   
  bNzero     <- nB!=0                    # T: if not all b's are zeros
  tb         <- array(0, c(nClass, nTest, nR))     
  overlapcstB  <- rep(0, nTest)
  overlapcstB1 <- rep(0, nTest)
  oltb         <- matrix(nClass, nTest)  # overlaped true b's
  for(i in 1:nR){
  tb[,,i]   <- cstB[[i]]
  overlapcstB  <- overlapcstB + (colSums(cstB[[i]])>0)*i
  overlapcstB1 <- overlapcstB1 + (colSums(cstB[[i]])>0)*1  
  }
  oltb <- apply(tb, c(1,2), sum)
  tb[tb==0]     <- nB + 1
  oltb[oltb==0] <- nB + 1
  if(!bNzero)  nR = 0;
  tr <- rep(1, nTest)
  
  if(max(overlapcstB1)>1){
  complexR <- TRUE
  oltb <- NULL
  }else{
  complexR <- FALSE
  tr       <-  overlapcstB 
  tr[tr==0] <- nR + 1
  }
  
  ta <- cstA                	
  ta[ta<0] <- nA + abs(ta[ta<0])
  # individual level data                 
  DAT <- matrix(0,N,nTest)
  for(j in 1:nTest) DAT[,j] <- rep(c(data[,j]), c(data[,nTest+1]))  
  patternInd <- rowSums(DAT*matrix(rep(2^((nTest-1):0), N), ncol=nTest, byrow=TRUE))+1
  
  # observed numbers in all patterns: begin
  C <- 2^nTest	
  allCell <- rep(0, C)                   # observed        
	for(i in 1:nCell){			    	
	index <- 1
	for(j in nTest:1){	
		index <- index + data[i, j]*2^(nTest-j)
		}		
	allCell[index] <- data[i, nTest+1]
	}		
	patterns <- rep("", C)
	for(c in 0:(C-1)){      # totally number of cells
	       	      num <- c
		      cell <- rep(0, nTest)      
		      for(j in nTest:1){        # P: the number of tests
		      cell[j] <- num - (num%/%2)*2          
		      num <- num%/%2                           
	      } 	      
	      patterns[c+1] <- paste(cell, collapse="")
	}
	# observed numbers in all patterns: end	
	
	patterns2 <- rep("", nrow(data))
	for(j in 1:nrow(data)){      # number of patterns in the data
	      patterns2[j] <- paste(data[j,1:nTest], collapse="")
	}	
	lcmr.obj <- list(N=N, nTest=nTest, test.names=test.names, patterns=patterns, patterns2=patterns2,   patternInd=patternInd, nCell=nCell, nClass=nClass, nA=nA, nB=nB, bNzero=bNzero, nR=nR, allCell=allCell, data=data, DAT=DAT, cstA=cstA, cstB=cstB, ta=ta, tb=tb, tr=tr, complexR=complexR, oltb=oltb)
	class(lcmr.obj) <- c("lcmr")
	lcmr.obj
}


#=================== print.lcmr ==========================
print.lcmr <- function(x, ...){
cat("Data set:\n")
print(x$data)
cat("Total number of observations: ", x$N, "\n")
cat("#classes:", x$nClass, "; #a's:", x$nA,  "#b's:", x$nB, "#r's:")
if(x$nB==0){
	cat("0\n")
}else{
cat(x$nR, "\n")
}
cat("Constraints on a's: \n")
print(x$cstA)

if(sum(x$nR)==0){
cat("There is no random effects.\n")
}else{
cat("Constraints on b's: \n")
print(x$cstB)
} 
}


