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


#------------------------------------------ mu --------------------------------------------------------------
mu <- function(x, start,end,...){
if(missing(end)) end<- x$iter
C      <- 2^x$nTest
burnin <- start-1
n.col  <- end-burnin
prev   <- x$prev[start:end, ]
a.tmp  <- cbind(x$a, 0)
a      <- a.tmp[start:end, ]
mu     <-  vector("list", n.col)

if(x$bNzero){
   b.tmp <- cbind(x$b, 0)
   b     <- b.tmp[start:end, ]
#  rArray <- x$rArray[,,start:end]   
  for(i in 1:n.col){
  muArray     <- array(0, c(x$N, x$nTest, x$nClass))
  for(j in 1:x$nClass){
  muArray[,,j]  <- pnorm(matrix(1, x$N, 1)%*%a[i, x$ta[j, ]] + x$rArray[,,start+i-1]%*%matrix(b[i, x$tb[j,,]], nrow=x$nR, byrow=TRUE));          # rArray: N*nR*newiter
}
  mu[[i]] <- muArray
}   # end of  for(i in 1:n.col)    
}else{
  for(i in 1:n.col){
  muArray     <- array(0, c(x$N, x$nTest, x$nClass))
  for(j in 1:x$nClass){
  muArray[,,j]  <- pnorm(matrix(1, x$N, 1)%*%a[i, x$ta[j, ]])
}
  mu[[i]] <- muArray
}   # end of  for(i in 1:n.col) 
} 
mu
}

