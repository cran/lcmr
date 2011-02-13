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

#corr <- function(x, ...)    UseMethod("corr")

#================== corr ========================
corr <- function(x, start, end, ...){

if(missing(end)){
end    <- x$iter
}

C      <- 2^x$nTest
nClass <- x$nClass
N      <- x$N  
nA     <- x$nA
nB     <- x$nB
nR     <- x$nR
tb     <- x$tb
ta     <- x$ta
burnin  <- start-1
n.col  <- end - burnin
prev   <- x$prev[start:end, ]
a.tmp  <- cbind(x$a, 0)
a      <- a.tmp[start:end, ]
bNzero <- x$bNzero

if(bNzero){
    b.tmp  <- cbind(x$b, 0)
    b      <- b.tmp[start:end, ]
}else{
    b <- matrix(0, nrow=n.col, ncol=1)
}

nCorr  <- (x$nTest-1)*x$nTest/2
corr   <- rep(0,  nrow=nCorr)
pCorr  <- matrix(0,  nrow=nCorr,  ncol=n.col)

# observed correlations
s     <- colSums(x$DAT)

k <- 0
for(i1 in 1:(x$nTest-1)){
    for(i2 in (i1+1):x$nTest){
        k <- k + 1
        corr[k] <- (sum(x$DAT[,i1]*x$DAT[,i2])*N - s[i1]*s[i2])/max(sqrt(s[i1]*s[i2]*(N - s[i1])*(N - s[i2])), 0.01)
    }
}

if(nR <=1){
tr     <- rep(1, x$nTest);
gq.out <- gauQuad()
gq     <- gq.out$gq
w      <- gq.out$w
for(i in 1:n.col){
# predicted (expected) correlations
s <- rep(0, x$nTest)
for(j in 1:x$nTest){
        for(n.i in 1:nClass){
        s[j] <- s[j] + prev[i, n.i]*pnorm(a[i,ta[n.i, j] ]/sqrt(1+b[i, tb[n.i, j, 1]]^2))
     }
 }


k <- 0
for(i1 in (1:(x$nTest-1))){
    for(i2 in ((i1+1):x$nTest)){
               k <- k+1
               if(tr[i1]==tr[i2]){
                       I <- 0
                       for(n.i in 1:nClass){
                        I <- I + prev[i, n.i]*sum(w*pnorm(a[i,ta[n.i, i1]]+ b[i, tb[n.i, i1, 1]]*gq)*pnorm(a[i, ta[n.i, i2]]+ b[i, tb[n.i, i2, 1]]*gq))
                       }
                       pCorr[k, i] <- (I-s[i1]*s[i2])/sqrt(s[i1]*s[i2]*(1-s[i1])*(1-s[i2]))
               }else{
                       I <- 0
                       for(n.i in 1:nClass){
                        I <- I + prev[i, n.i]*pnorm(a[i,ta[n.i, i1]]/sqrt(1+b[i, tb[n.i, i1, 1]]^2))*pnorm(a[i,ta[n.i, i2]]/sqrt(1+b[i, tb[n.i, i2, 1]]^2))
                       }
                       pCorr[k, i] <- (I-s[i1]*s[i2])/sqrt(s[i1]*s[i2]*(1-s[i1])*(1-s[i2]))
                       }
               }
}

}

}else{
nrnd <- 1000;
nallrnd <- nrnd*nR;
for(i in 1:n.col){
# predicted (expected) correlations
s <- rep(0, x$nTest)
for(j in 1:x$nTest){
for(n.i in 1:nClass){
    s[j] <- s[j] + prev[i, n.i]*pnorm(a[i,ta[n.i, j] ]/sqrt(1+sum(b[i, tb[n.i, j, ]]^2)))
     }
}       # end of for(j in 1:x$nTest)


rndMat <- matrix(rnorm(nallrnd), nrnd, nR)        # generating random numbers. Do I need to generate them for each iteration?
pd     <- matrix(0,  nrow=nrnd, ncol=nClass)
for(n.i in 1:nClass){
    sumabrs1 <- a[i, ta[n.i, i1]]
    sumabrs2 <- a[i, ta[n.i, i2]]
    for(n.r in 1:nR){
      sumabrs1 <- sumabrs1 + b[i, tb[n.i, i1,n.r]]*rndMat[,n.r]
      sumabrs2 <- sumabrs2 + b[i, tb[n.i, i2,n.r]]*rndMat[,n.r]
    }
    pd[, n.i] <- prev[i, n.i]*(pnorm(sumabrs1)*pnorm(sumabrs2))
}

k <- 0
for(i1 in (1:(x$nTest-1))){
   for(i2 in ((i1+1):x$nTest)){
   k <- k+1
   I <- 0
   for(n.i in 1:nClass){
        sumabrs1 <- a[i, ta[n.i, i1]]
        sumabrs2 <- a[i, ta[n.i, i2]]
        for(n.r in 1:nR){
            sumabrs1 <- sumabrs1 + b[i, tb[n.i, i1,n.r]]*rndMat[,n.r]
            sumabrs2 <- sumabrs2 + b[i, tb[n.i, i2,n.r]]*rndMat[,n.r]
        }
        I <- I + prev[i, n.i]*mean(pnorm(sumabrs1)*pnorm(sumabrs2))
    }
    pCorr[k, i] <- (I-s[i1]*s[i2])/sqrt(s[i1]*s[i2]*(1-s[i1])*(1-s[i2]))
    }
}
}    # end of for(i in 1:n.col)
}    # end of if-else
list(corr=corr, pCorr=pCorr)
}


