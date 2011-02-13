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

#------------------------------------------ probT1 --------------------------------------------------------------
probT1 <- function(x, start, end, ...){

if(missing(end)) end <- x$iter 
burnin <- start-1

n.col  <- end-burnin
prev   <- x$prev[start:end, ]
a.tmp  <- cbind(x$a, 0)
a      <- a.tmp[start:end, ]
if(x$bNzero){
b.tmp  <- cbind(x$b, 0)
b      <- b.tmp[start:end, ]
}else b <- matrix(0, nrow=n.col, ncol=1)

probT1  <- array(0, c(x$nClass, x$nTest, n.col))
for(i in 1:n.col){
    for(n.i in 1:x$nClass){
    sumbs <- 0
    if(x$nR>0){
    for(n.r in 1:x$nR){
     sumbs <- sumbs + b[i, x$tb[n.i,,n.r]]^2
    }
    }
probT1[n.i,,i]  <- pnorm(a[i, x$ta[n.i, ]]/sqrt(1+ sumbs))
}
}
probT1
}
