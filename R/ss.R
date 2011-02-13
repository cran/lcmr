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


#==================ss =======================================
ss <- function(x, start, end, D=1, probs= c(0.025, 0.50, 0.975), ...){     # sensitivities and  specificities
if(missing(end)) end <- x$iter

burnin  <- start-1

if(is.null(x$probT1)){
  x$probT1 <- probT1(x,  start, end)
}
nClass     <- x$nClass
nTest      <- x$nTest
test.names <- x$test.names
n.iter     <- end-burnin              # number of iterations for the calculation
sen        <- matrix(0, n.iter, nTest)    # sensitivities
p.i        <- array(0, c(n.iter, nTest, nClass))

for(i in D){
    p.i[,,i] <- matrix(rep(x$prev[start:end, i], nTest), ncol=nTest, byrow=FALSE)
    sen <- sen + t(x$probT1[i, ,])*p.i[,,i]
}

sen.tmp <- sen/apply(p.i, c(1,2), sum)
lenProb <- length(probs)
mat1 <- matrix(0, nTest, (1+lenProb))
colnames(mat1) <- c(paste(probs*100, "%"), "Mean")
# rownames(mat1) <- paste("Test", as.character(1:nTest), sep="")
rownames(mat1) <- test.names

for(j in 1:nTest){
                mat1[j, 1:lenProb] <- quantile(sen.tmp[, j], probs=probs, na.rm = TRUE)
                mat1[j, 1+lenProb] <- mean(sen.tmp[, j],  na.rm = TRUE)
    }



p.i <- array(0, c(n.iter, nTest, nClass))
spec <- matrix(0, n.iter, nTest)

D1 <- rep(0, nClass)
D1[D] <- D
tmp <- c(1:nClass)
for(i in tmp[!tmp==D1]){
p.i[,,i] <- matrix(rep(x$prev[start:end, i], nTest), ncol=nTest, byrow=FALSE)
spec <- spec + (1-t(x$probT1[i,,]))* p.i[,,i]
}

spec.tmp <-spec/apply(p.i, c(1,2), sum)

mat2 <- matrix(0, nTest, (1+lenProb))
colnames(mat2) <- c(paste(probs*100, "%"), "Mean")

# rownames(mat2) <- paste("Test", as.character(1:nTest), sep="")
rownames(mat2) <- test.names
for(j in 1:nTest){
                mat2[j, 1:lenProb] <- quantile(spec.tmp[, j], probs=probs, na.rm = TRUE)
                mat2[j, 1+lenProb] <- mean(spec.tmp[, j],  na.rm = TRUE)
    }

list(sensitivity=mat1, specificity=mat2)
}

