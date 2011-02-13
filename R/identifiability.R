                                

# check identifiability by calculating Jacobian
identifiability <- function(lcmrObj){
	if(class(lcmrObj)!="lcmr") stop("The input argument is an object of the class lcmr!")
	    
#------ variables coming form lcmrObj-----------
K        <- lcmrObj$nCell
nTest    <- lcmrObj$nTest
nPattern <- 2^nTest
ndf      <- nPattern -1
N        <- lcmrObj$N
allCell  <- lcmrObj$allCell
nA       <- lcmrObj$nA
nB       <- lcmrObj$nB
nR       <- lcmrObj$nR
bNzero   <- (nB!=0)                         # T: if not all b's are zeros
nClass   <- lcmrObj$nClass
nPar     <- nClass -1 + nA + nB # number of parameters
tt       <- lcmrObj$data
tb       <- lcmrObj$tb 
ta       <- lcmrObj$ta 

a     <- c(runif(nA), 0)
b     <- c(abs(runif(nB)), 0)
prev  <- runif(nClass)
prev  <- prev/sum(prev)

nAs <- matrix(0, nrow=nClass, ncol=nA)          
for(i in 1:nClass){
    for(j in 1:nA){        
    nAs[i, j] <- sum((ta[i, ]== j))   
    }
}

# all possible patterns
patterns <- matrix(0, nPattern, nTest)	
for(c in 0:(nPattern-1)){      # totally number of cells
	num <- c		      
	for(j in nTest:1){        # P: the number of tests
		patterns[c+1, j] <- num - (num%/%2)*2          
		num <- num%/%2                           
	}   	      
}             
                 
#----generate the random effects----------------
M <- 1000
z <- matrix(rnorm(M*nR), M, nR) 
sumabrsArray <- array(0, c(M,  nClass, nTest));		
if(nR>0){
		for(i in 1:nTest){			
     		  for(n.i in 1:nClass){
  		      sumabrsArray[, n.i, i] <- a[ta[n.i, i]] + z%*%matrix(b[tb[n.i, i, ]], ncol=1);   
		        } # end of for(n.i in 1:nClass)
		} # end of for(i in 1:nTest)
		}else{
    for(i in 1:nTest){			
     		  for(n.i in 1:nClass){
  		      sumabrsArray[, n.i, i] <- a[ta[n.i, i]];   
		        } # end of for(n.i in 1:nClass)
		} # end of for(i in 1:nTest)
		}
		

muabrArray <- pnorm(sumabrsArray)	
gArray  <- array(0, c( M, nTest, nClass, nPattern))
for(k in 1:nClass){
   for(j in 1:M){
    for(i in 1:nPattern){
     gArray[j, , k, i] <- muabrArray[j, k, ]^patterns[i,]*(1-muabrArray[j, k, ])^(1-patterns[i,])    
    }
   }
}


#---- partial derivatives ----------------	 	 
partialMat <- matrix(0, nPattern, nPar) 
for(i in 1:nPattern){
  for(k in 1:nClass-1){
  partialMat[i, k] <- (sum(apply(gArray[, ,k, i], 1, prod)) - sum(apply(gArray[,,nClass, i], 1, prod)))/M
  }
  tmpDeva <- rep(0, nA+1) 
  tmpDevb <- rep(0, nB+1) 
  for(k in 1:nClass){
     for(l in 1:nTest){      
       tmpDeva[ta[k, l]] <- tmpDeva[ta[k, l]] + prev[k]*(2*patterns[i, l]-1)*sum(dnorm(sumabrsArray[, k, l])*muabrArray[,  k, l]* apply(gArray[, l!=c(1:nTest), k, i], 1, prod))       
       if(nB>0){
       for(q in 1:nR){       
       tmpDevb[tb[k, l, q]] <- tmpDevb[tb[k, l, q]] + prev[k]*(2*patterns[i, l]-1)*sum(z[,q]*dnorm(sumabrsArray[, k, l])*muabrArray[,  k, l]* gArray[, l!=c(1:nTest), k, i])       
       }
       }
     }     
  }
  partialMat[i, nClass:(nClass+nA-1)] <- tmpDeva[1:nA]/M
  if(nB>0) partialMat[i, (nClass+nA):nPar] <- tmpDevb[1:nB]/M 
}
	qrrank = qr(partialMat)$rank		
	identifiable = 0           # not identifiable
	
	# necessary condition: the number of parameters <= number of degrees of freedom	
	# sufficient condition: rank of the Jacobian ==  the number of parameters
		
	if(nPar<=ndf){		
		if(qrrank==nPar){	
		identifiable <- 1          # identifiable
		}else{		
		identifiable <- 2          # might not identifiable
		}
	}	
			
	cat("Number of parameters: ", nPar, "\nDegrees of freedom:", ndf, "\n")	
	cat("Rank of the Jacobian: ", qrrank, "\n")	
	if(identifiable==0) cat("The model is not identifiable. \n")
	if(identifiable==1) cat("The model is identifiable. \n")
	if(identifiable==2) cat("The sufficient condition of identifiability is not satisfied. The model might not be identifiable. \n")	
	jacob <- list(nPar=nPar, ndf=ndf, rank=qrrank, identifiable=identifiable)
}

