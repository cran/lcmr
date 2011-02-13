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
function(object, ...)    UseMethod("logLik")
#------------------------------------------ logLik  --------------------------------------------------------------
logLik.gibbs <- function(object, start,end,...){

if(missing(end))  end <- object$iter

n.col  <- end-start+1
prev   <- object$prev[start:end, ]
a.tmp  <- cbind(object$a, 0)
a      <- a.tmp[start:end, ]
logLik <- rep(0, n.col)

uniquePn <- unique(object$patternInd)
nUniqPn  <- length(uniquePn)


if(object$bNzero){
b.tmp <- cbind(object$b, 0)
b     <- b.tmp[start:end, ] 

pd.i.sum<-rep(0,nUniqPn)
create.seq<-function(x) {seq(1,20)}
nRgrid<-expand.grid(lapply(seq(1,object$nR),create.seq))

library(statmod)
out <- gauss.quad(20,"hermite")
gauss.quad.nodes <- sqrt(2)*as.vector(out$nodes)   
gauss.quad.wts <- as.vector(out$weights)/sqrt(pi) 

wt.matrix<-matrix(0,dim(nRgrid)[1],object$nR) 
for (i in 1:object$nR) wt.matrix[,i]<-gauss.quad.wts[nRgrid[,i]]
gwts<-apply(wt.matrix,1,prod)

node.matrix<-matrix(0,dim(nRgrid)[1],object$nR) 
for (i in 1:object$nR) node.matrix[,i]<-gauss.quad.nodes[nRgrid[,i]]
prob<-matrix(0,nUniqPn,n.col)
for(i in 1:n.col){
    pd.i    <- matrix(0,nUniqPn,dim(nRgrid)[1])
    for(ii in 1:dim(nRgrid)[1]){
    for(j in 1:object$nClass){
  
    muMat1<-matrix(1, nUniqPn, 1)%*%a[i, object$ta[j, ]]
    muMat2<-matrix(rep(node.matrix[ii,],nUniqPn),ncol=object$nR,byrow=T) %*% matrix(b[i, object$tb[j,,]], nrow=object$nR, byrow=TRUE)
    
      muMat     <- pnorm(muMat1+muMat2);          # rArray: N*nR*newiter
      pd.i[, ii] <- pd.i[,ii]+ prev[i, j]*apply(muMat*object$data[,1:object$nTest] + (1-muMat)*(1-object$data[,1:object$nTest]), 1, prod);
    }    
    pd.i[,ii]<-pd.i[,ii]*gwts[ii]
}  

prob[,i]<-apply(pd.i,1,sum)
logLik[i] <- sum(object$data[,object$nTest+1]*log(prob[,i])) 
}   # end of  for(i in 1:n.col)
dhat.test2<--2*sum(object$data[,object$nTest+1]*log(rowMeans(prob)))


}else{
amean    <- c(colMeans(object$a[start:end, 1:object$nA]), 0); 
prevmean <- colMeans(object$prev[start:end, ]); 

logLik <- object$logLik[start:end]

nrowData <- nrow(object$data);	
pd.i    <- matrix(0,  nrowData,  object$nClass)     
for(j in 1:object$nClass){
    pnormabr  <- pnorm(matrix(1, nrowData, 1)%*%matrix(amean[object$ta[j, ]],nrow=1)); 
    pd.i[, j] <- prevmean[j]*apply(pnormabr*object$data[,1:object$nTest]+(1-pnormabr)*(1-object$data[,1:object$nTest]), 1, prod); 
    }  
prob<-apply(pd.i,1,sum)
dhat.test2=-2*sum(object$data[,object$nTest+1]*log(rowSums(pd.i)))
}

dbar<--2*mean(logLik)
pd.test2<-dbar-dhat.test2
dic.test2<-dhat.test2+2*pd.test2

sumy.logLik <- quantile(logLik, c(0.025, 0.5, 0.975))

   list(sumy.logLik = sumy.logLik, logLik=logLik,dbar = dbar, dhat = dhat.test2, 
        pd = pd.test2, dic = dic.test2,prob=prob)

}


