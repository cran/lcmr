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

#------------------------------------------ probD  --------------------------------------------------------------
#P(D|profile)
probD <- function(x,  start, end, probs=c(0.025, 0.50, 0.975), ...) {
if(missing(end)) end <- x$iter

#ptm = proc.time()
#cat("Calculating Pr(latent class|combination of observed variables)...\n")

#burnin  <- start - 1
#n.col   <- end-burnin
#prev    <- x$prev[start:end, ]

#if(!is.null(x$mu)){
#muList  <- x$mu;
#}else{
#muList  <- mu(x, burnin);
## x$mu    <- muList;
#}
# muArray, N*nTest*nClass*n.col
#cellMat <- x$data[,1:x$nTest] 
#n.dcell <- nrow(cellMat)
#probD   <- array(0,  c(n.dcell,x$nClass, n.col))

#for(i in 1:n.col){
#randProd <- matrix(1, x$N, x$nClass)
#for(n.i in 1:x$nClass){	 
#      if(sum(x$whichR[n.i,])>0){
#       randProd[, n.i]  <- apply(dnorm(as.matrix(x$rArray[,(1:x$nR)[x$whichR[n.i,]],burnin+1])), 1, prod)
#	  }
#}
#cSrProd <- colSums(randProd)

#pd      <- matrix(0,  nrow=x$N, ncol=x$nClass)
#muArray <-  muList[[i]]
#for(j in 1:n.dcell){      # totally number of cells
#	 curCell <- matrix(rep(as.numeric(cellMat[j,]), x$N), nrow=x$N, byrow=TRUE)
#	 for(n.i in 1:x$nClass){	 
#     pd[, n.i] <- prev[i, n.i]*apply(muArray[,,n.i]*curCell + (1-muArray[,,n.i])*(1-curCell), 1, prod)*randProd[,n.i]/cSrProd[n.i]
#     }
#     probD[j, ,i] <- colSums(pd)/sum(pd)
#    } # end of for(j in 1:n.dcell)
#} # end of  for(i in 1:n.col)
#print(proc.time() -ptm)

probD <- x$probD[,,start:end]
str1  <- "list("
for(i in 1:(x$nClass-1)) str1 <- paste(str1, "Class", i, "=NULL,", sep="")
str1  <- paste(str1, "Class", x$nClass, "=NULL)", sep="")
summary.probD <- eval(parse(text=str1)) 
quantArray <- apply(probD, c(1, 2), quantile,  probs=probs, na.rm = TRUE)
meanMat    <- apply(probD, c(1, 2), mean, na.rm = TRUE)
for(i in 1:x$nClass){      
        mat <- t(rbind(quantArray[,,i], mean=meanMat[,i]))            
        rownames(mat)    <- paste("P(Class ", i, "|", x$patterns2,")", sep="")                        
        summary.probD[i] <- list(mat)
      }      
summary.probD	               
}

