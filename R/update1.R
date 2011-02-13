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

#============== update.gibbs ======================================
#------ main function of the gibbs sampler-----------------
update1  <- function(object, newiter, ...){
# object is an object of the class gibbs
ptm <- proc.time()
prior.sdA <- object$prior.sdA
prior.mA  <- object$prior.mA
sdA       <- object$sdA
mA        <- object$mA
priorPrev <- object$priorPrev
priorA    <- object$priorA
priorB    <- object$priorB
bNzero    <- object$bNzero                    # T: if not all b's are zeros
nTest     <- object$nTest
N         <- object$N
allCell   <- object$allCell
nA        <- object$nA
nB        <- object$nB
nR        <- object$nR
nClass    <- object$nClass
tb        <- object$tb;
ta        <- object$ta;
data      <- object$data
DAT       <- object$DAT      # individual level data
last      <- object$last
init      <- object$init
tr        <- object$tr
complexR  <- object$complexR

C         <- 2^nTest
iter      <- last + newiter
#-------------------------------------------------------------------------------------
nAs <- matrix(0, nrow=nClass, ncol=nA)
for(i in 1:nClass){
     for(j in 1:nA){
       nAs[i, j] <- sum((ta[i, ]== j))
     }
}

# defining variables for the gibbs sampler
newiter  <- newiter + 1
prev     <- matrix(0, newiter, nClass)
negA     <- object$cstA[object$cstA<0]
a        <- matrix(0, newiter, (nA + length(unique(negA))))
if (length(unique(negA))>0){
        a[,(nA+1):(nA+length(unique(negA)))] <- prior.mA[(nA+1):(nA+length(unique(negA)))]
}
b        <- matrix(0, newiter, (nB + 1))
N.d      <- rep(0, nClass)

# probD
uniqPns <- unique(object$patternInd)
probD   <- array(0,  c(length(uniqPns),nClass, iter))
if(last>0){
   probD[,,1:last] <- object$probD
}

# for the purpose of posterior predictive check
texp    <- matrix(0, N, nTest)
predmat <- matrix(0, newiter, 2^nTest)

  nCorr   <- nTest*(nTest-1)/2
  agree   <- array(0,  c(nCorr, nClass, iter))   #  observed agreement between two tests
  pAgree  <- array(0,  c(nCorr, nClass, iter))   # predicted agreement between two tests
  if(last>0){
  agree[,,1:last]  = object$agree$agree;
  pAgree[,,1:last] = object$agree$pAgree;
  }

logLik  <- rep(0, newiter);
pd.i    <- matrix(1,  N,  nClass)

# starting values
if(last>0){
a[1, 1:nA]  <- object$a[last,1:nA]
if(bNzero)  b[1, 1:nB]  <- object$b[last, ]
prev[1, ]   <- object$prev[last, ]
}else{
  if(bNzero){
    b[1, 1:nB] <- init$b
   }
a[1, 1:nA]  <- init$a
prev[1, ]   <- init$prev
}

# r <- matrix(rep(rnorm(N), nR), nrow=N)       # initial random effects
r <- cbind(object$lastR, 0)                    # last random effects from last run.


# gibbs sampler iterations - order follows Qu & Hadgu, proceedings paper
if(last>0){
y <- object$lastY
}

for(i in 2:newiter){
# Step 1
for(j in 1:nClass){
    A <- matrix(rep(a[i-1, ta[j,]], N), nrow=N, byrow=TRUE)
    pnormA <- pnorm(A)
    pd.i[,j]  <- prev[i-1, j]*apply((pnormA*DAT+(1-pnormA)*(1-DAT)), 1, prod)
}
# logLik
logLik[i]  <- sum(log(rowSums(pd.i)))

pd.i <- pd.i/matrix(rep(rowSums(pd.i), nClass), ncol=nClass)

# probD
for(j in 1:length(uniqPns)){
probD[j,,i+last-1] <- pd.i[which(object$patternInd==uniqPns[j])[1],]
}


d <- apply(pd.i, 1, rmulti, size=nClass)
for(j in 1:nClass) N.d[j] <- sum(d==j)
D <- matrix(rep(d, nTest), ncol=nTest)

if(last==0 & i==2){
    # each row are propabilities, which are internally normalized to sum 1.
    # d <- apply(pd.i, 1, rmulti, size=nClass)
        y <- 2*DAT-1
        if(length(priorA$u)==nA | length(priorA$l)==nA){
            indx.a <- (1:nA)[priorA$u < 1000 & priorA$l > -1000]
        for(a.j in indx.a){
        mid.value <- (priorA$u[a.j]-priorA$l[a.j])/2
        pos <- (1:(nClass*nTest))[c(ta==a.j)]
        for(p.j in pos){
        col.pos <-  ceiling(p.j/nClass)
        row.pos <- p.j - (col.pos-1)*nClass
        y[d==row.pos, col.pos] <-  mid.value
        }
        }
        }
}

# Step 2  for r # no this step
# Step 3 for b, # no this step

# Step 4
if(nClass==2){
    prev[i, 1] <- rbeta(1, priorPrev[1] + N.d[1],  priorPrev[2] + N.d[2])
    prev[i, nClass] <- 1 - sum(prev[i,1])
}else{
    prev[i,] <- rdirichlet(1, c(priorPrev + N.d))
}


# Step 5
m.y <- matrix(0, nrow=N, ncol=nTest)
for(n.i in 1:nClass){
    for(j in 1:nTest){
    m.y[d==n.i, j]  <-  a[i-1, ta[n.i, j]]
    }
}

y[DAT==1 & m.y > 8]  <- 8
y[DAT==0 & m.y < -8] <- -8
y[DAT==1 & m.y < 8]  <- truncnorm(sum(DAT==1 & m.y < 8),  m.y[DAT==1 & m.y <  8], 1, 0, 8)
y[DAT==0 & m.y > -8] <- truncnorm(sum(DAT==0 & m.y > -8), m.y[DAT==0 & m.y > -8], 1, -8, 0)


#Step 6
m.a  <- rep(0, nA+1)
v.a  <- rep(0, nA+1)
v.a[1:nA] <-  1/((1/sdA[1:nA]^2)+ N.d %*% nAs)
m.a[1:nA] <- (mA[1:nA]/(sdA[1:nA]^2))

for(n.i in 1:nClass){
    for(j in 1:nTest){
    m.a[ta[n.i, j]] <- m.a[ta[n.i, j]] +  sum(y[d==n.i, j])
    }
}

m.a[1:nA]   <- v.a[1:nA]*m.a[1:nA]
a[i,1:nA]   <- truncnorm(n=nA, m.a[1:nA], v.a[1:nA], priorA$l, priorA$u)

a[i,a[i,]==Inf] <-0.5*(priorA$l[a[i,]==Inf]+priorA$u[a[i,]==Inf])
a[i,a[i,]==-Inf]<-0.5*(priorA$l[a[i,]==-Inf]+priorA$u[a[i,]==-Inf])

# predicted observations
for(j in 1:nTest){
            prob <- rep(0, N)
            for(n.i in 1:nClass){
            prob[d==n.i] <- pnorm(a[i, ta[n.i, j]])
            }
            texp[, j] <- rbinom(N, 1, prob)
    }
cellposVec <- texp%*%(2^(seq(nTest-1,0)))+1
for(n.c in 1:C){
    predmat[i, n.c] <- sum(cellposVec==n.c)
}

k <- 0
for(i1 in 1:(nTest-1)){
            for(i2 in (i1+1):nTest){
            k <- k + 1
            for(n.i in 1:nClass){
            pAgree[k, n.i, i+last-1]  <- sum(texp[d==n.i, i1]== texp[d==n.i, i2])/N.d[n.i]
            agree [k, n.i, i+last-1]  <- sum(DAT[d==n.i, i1] == DAT[d==n.i, i2])/N.d[n.i]
            }
        }
    }

} #end of for (i in 2:newiter)



if (length(unique(negA))>0) {
   object$a  <-  rbind(object$a, a[2:newiter, 1:(nA+length(unique(negA)))])
}else{
   object$a  <-  rbind(object$a, a[2:newiter, 1:nA])
}

if(bNzero){
   if(nB==1){
    object$b <- rbind(object$b, matrix(b[2:newiter,1], ncol=1))
    }else{
    object$b <- rbind(object$b, b[2:newiter, 1:nB])
    }
}else{
    object$b <- NULL
}

object$logLik         <-  c(object$logLik, logLik[2:newiter])
object$prev           <-  rbind(object$prev, prev[2:newiter, ])
if(nR>0) object$lastR <- r[, 1:nR]
object$lastY          <- y
object$probD          <- probD
object$agree$agree    <-  agree
object$agree$pAgree   <-  pAgree
object$predmat        <- rbind(object$predmat, predmat[2:newiter,])
object$iter           <- object$iter + newiter - 1
object$last           <- object$iter
mininutes <- (proc.time() - ptm)[3]/60
cat("Iterations of the Gibbs sampler:", newiter-1, "\nTime elapsed:", mininutes, "mininutes. \n")
return(object)
}
