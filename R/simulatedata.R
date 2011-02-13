# simulate some data sets.	
simuData <- function(N,  prev, trueP){
trueP  <- t(trueP)
P	   <- nrow(trueP)
nClass <- ncol(trueP)
simu.data <- matrix(0, N, P)
d         <- rep(0, N)
N.d       <- rep(0, nClass)

for(j in 1:N){
        d[j] <- rmulti(prev, nClass)
	}

for(i in 1:nClass){
	N.d[i] <- sum(d==i)
	}

for(j in 1:nClass){
	for(i in 1:P){
		simu.data[d==j, i] <- rbinom(N.d[j], 1, trueP[i, j])				
	}
}

cell <- 1
for(i in P:1){
	cell <- cell + simu.data[, i]*2^(i-1)
}

C  <- sort(unique(cell))
len.C <- length(C)
dd <- matrix(0, len.C, P+1)   # data set
colnames(dd) <- c(paste("Test", c(1:P), sep=""), "count")

for(i in 1:len.C){
    dd[i, (P+1)] <- sum(cell==C[i])
    if(sum(cell==C[i])>1)    dd[i, 1:P] <- simu.data[cell==C[i], ][1, ]
    else{
      dd[i, 1:P] <- simu.data[cell==C[i], ]
    } 
}
                                                       
p <- matrix(0, P, nClass)
for(j in 1:nClass){
	for(i in 1:P){
		p[i, j] <- sum(simu.data[d==j, i])/N.d[j]																													
		
	}
}

return(list(dat=dd, prev=N.d/N, p=p))
}





compactMatrix <- function(M, cols){                                   # M is a matrix , cols: columns selected
n.col <- ncol(M)
n.row <- nrow(M)
n     <- length(cols)+1     # number of columns
newM  <- M[, c(cols, n.col)]
if(is.null(colnames(newM)))
colnames(newM) <- c(paste("Test", cols, sep=""), "count")

for(i in 1:(n.row-1)){
        if(newM[i, n] > 0){
	for(j in (i+1):n.row){
	        if(newM[j, n] > 0){	     	        	        
		if(sum(M[j, cols]==M[i, cols])==(n-1)){
		newM[i, n] <- newM[i, n] + newM[j, n]
		newM[j, n] <- -100
		}
		}
	}
	}
}
return(newM[newM[, n]>0, ])
}




