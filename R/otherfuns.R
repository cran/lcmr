#--------------------------------------------------------------------------------------

parRC <- function(Q){
num.r <- ceiling(sqrt(Q))
num.c <- num.r
num <- num.r * num.c
min1 <- max(1, num.r-1)
max1 <- min(Q, num.r+1)
min2 <- max(1, num.c-1)
max2 <- min(Q, num.c+1)
for(i in min1:max1){
    for(j in min2:max2){
      new.num <- i*j
if((new.num>=Q) & (new.num < num)){
num <- new.num
num.r <- i
num.c <- j
}
    }
}
    return(c(min(num.r, num.c), max(num.r, num.c)))
}


#--------------------------------------------------------------------------------------

gauQuad <- function(){
    # default number of the gaussian quadrature values is 20
    numGQ <- 20
    gq <- c(-7.619049, -6.51059, -5.578739,-4.734581,-3.943967,-3.189015,-2.458664,-1.745247,-1.042945,-0.3469642,0.3469642,1.042945,1.745247,2.458664,3.189015,3.943967,4.734581,5.578739,6.51059,7.619049)
    w  <- c(1.257801e-13,2.482062e-10,6.12749e-08,4.402121e-06,0.0001288263,0.001830103,0.01399784,0.06150637,0.1617393,0.2607930,0.2607930,0.1617393,0.06150637,0.01399784,0.001830103,0.0001288263,4.402121e-06,6.12749e-08,2.482062e-10,1.257801e-13)
    list(numGQ=numGQ, gq=gq, w=w)
}

#------------------------------------------ rmulti  --------------------------------------------------------------

rmulti <- function(x, size) {
    return(sum(c(1:size)*rmultinom(1,1,x)))
}


#------------------------------------------ truncnorm  --------------------------------------------------------------
# function to sample from a truncated normal
truncnorm <- function(n=1,m,v,l,u) {
l1 <- rep(0, n)
u1 <- rep(1, n)
trunc.l <- l!=-1000        # index for the lower limit for the truncted parameters
trunc.u <- l!=1000         # index for the upper limit for the truncted parameters
l1[trunc.l]  <- pnorm((l-m)/sqrt(v))[trunc.l]
u1[trunc.u]  <- pnorm((u-m)/sqrt(v))[trunc.u]
           x <- runif(n,l1,u1)
           y <- qnorm(x)*sqrt(v)+m
           return(y)
    }


#============ rdirichlet ===================================
## pick n random deviates from the Dirichlet function with shape
## parameters a
rdirichlet<-function(n, a){
        l<-length(a);
        x<-matrix(rgamma(l*n,a),ncol=l,byrow=TRUE);
        sm<-x%*%rep(1,l);
        x/as.vector(sm);
    }

