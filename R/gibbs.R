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
#  http://www.r-project.org/Licenses


gs <- 
## Gibbs sampler
   function(lcmr.obj,
            iter,
			init,  
			prior.sdA, 
			prior.mA, 
			priorPrev,  
			priorA=list(l=-1000, u=1000), 
			priorB=list(l=0, u=5),...)
{ 
    call <- match.call();
   #--------------------check the input arguments----------------
    if(lcmr.obj$nClass!=length(init$prev)) 
           stop("The number of initial values for the prevalences is not correct!")
    if(lcmr.obj$nA!=length(init$a))
           stop("The number of initial values for the parameter a's is not correct!")
    if(lcmr.obj$nB!=length(init$b) & lcmr.obj$nB!=0){           
    warning("The number of initial values for the parameter b's is not correct. All b's are initialized 1!")                            
    init$b <- rep(1, lcmr.obj$nB)
    }   
    if(iter<0)  stop("The number of iterations cannot be negative!")
    
    # checking input parameters
    nClass <- lcmr.obj$nClass
    nA     <- lcmr.obj$nA
    nB     <- lcmr.obj$nB   
    if(missing(prior.sdA)){
#    warning("prior.sdA is missing! Using the default values (all ones)!\n")
    prior.sdA <- rep(1, nA)
    }else{
    if ((length(prior.sdA)!=length(unique(c(lcmr.obj$cstA)))) || (sum(prior.sdA<0)>0)) {
#    warning("The length of prior.sdA is not the same as the number of a's, or prior.sdA is negative! Using the default values!\n")
    prior.sdA <- rep(1, nA)
    }  
    }
    sdA <- c(prior.sdA, 1)
    
    if(missing(prior.mA)){
 #   warning("prior.mA is missing! Using the default values (all zeros)!\n")
    prior.mA <- rep(0, nA)
    }else{
    if (length(prior.mA)!=nA) {
 #   warning("The length of prior.mA is not the same as the number of a's! Using the default values!\n")
    prior.mA <- rep(0, nA) 
    }
    }
    mA <- c(prior.mA, 0)
    if(missing(priorPrev)){
     priorPrev <- rep(1, nClass)
     }else{
    if(length(priorPrev)!=nClass){
    warning("The length of priorPrev is not the same as the number of classes! Using the default values!\n")
    priorPrev <- rep(1, nClass)
    } 
    }

#--------------------initialize the Gibbs Sampler----------------
    pd.i <- matrix(0.5, lcmr.obj$N,  lcmr.obj$nClass)            
    d    <- apply(pd.i, 1, rmulti, size=lcmr.obj$nClass)        
    r    <- matrix(rnorm(lcmr.obj$N*lcmr.obj$nR),  nrow=lcmr.obj$N)                       # initial random effects                         
    gs.obj <- list(call=call, N=lcmr.obj$N, nTest = lcmr.obj$nTest, test.names = lcmr.obj$test.names,
    patterns = lcmr.obj$patterns,    patterns2 = lcmr.obj$patterns2,
	patternInd = lcmr.obj$patternInd, nCell = lcmr.obj$nCell, nClass=lcmr.obj$nClass,
    nA=lcmr.obj$nA, nB=lcmr.obj$nB, bNzero=lcmr.obj$nB, nR=lcmr.obj$nR, allCell=lcmr.obj$allCell,
    data=lcmr.obj$data, DAT=lcmr.obj$DAT, cstA =lcmr.obj$cstA, cstB=lcmr.obj$cstB,
    ta=lcmr.obj$ta, tb=lcmr.obj$tb, tr=lcmr.obj$tr, complexR=lcmr.obj$complexR,
    oltb = lcmr.obj$oltb,  iter=0,  init=init,
      prior.sdA=prior.sdA, sdA = sdA, prior.mA=prior.mA, mA=mA, 
      priorPrev=priorPrev, priorA=priorA, priorB=priorB, prev=NULL, a=NULL, b=NULL, lastR=r, lastD=d,last=0)        
    if(gs.obj$nR>0){
	tmpLogic       <- gs.obj$tb>0 & gs.obj$tb<=gs.obj$nB; 	
    gs.obj$whichR  <- apply(tmpLogic, c(1,3), sum)>0;
	}else{
	gs.obj$whichR  <- NULL;
	}                                                                                                                        
    class(gs.obj)  <- c("gibbs")

    
    if(iter>0){ 
    gs.obj  <- update(gs.obj, iter)                  # main part of the result of gibbs sampler                                                    
    } # end of if(iter>0)
cat("The total number of iteration is:", gs.obj$iter, ". Please use function 'summary' to summarize the results.\n")
return(gs.obj)
}

#=================== print.gibbs =============================
print.gibbs <- function(x,...){     
    burnin <- ceiling(0.1*x$iter)
    prev   <- rep(0, x$nClass)        
    names(prev) <- paste("Class", as.character(1:x$nClass))           
    for(j in 1:x$nClass) prev[j] <- mean(x$prev[burnin:x$iter, j],  na.rm = TRUE)           
    
    a <- rep(0, x$nA) 
    names(a) <- paste("a", as.character(1:x$nA), sep="")
    
    for(j in 1:x$nA)  a[j] <- mean(x$a[burnin:x$iter, j],  na.rm = TRUE)        
    cat("Iteration: ", burnin, "--", x$iter, "\n")    
    cat("Prevalence:\n")
    print(round(prev, 3))
    cat("\n Parameter a's:\n")
    print(round(a, 2))
    if(x$bNzero){
        b <- rep(0,x$nB)                  # judge whether b exists    
        names(b) <- paste("b", as.character(1:x$nB), sep="")
        for(j in 1:x$nB)  b[j] <- mean(x$b[burnin:x$iter, j],  na.rm = TRUE)            
        cat("\n Parameter b's:\n")
        print(round(b, 2))  
    }   
}

#====================== update.gibbs======================================
update.gibbs  <- function(object, newiter,  ...){
# object is an object of the class gibbs
if(object$nR ==0){
   object <- update1(object, newiter);
}else{
   if((object$nR ==1)|!object$complexR){
   object <- update2(object, newiter);
   }else{
   object <- update3(object, newiter);
   }
}
}


#======================summary.gibbs======================================
summary.gibbs <- function(object, start, end, D, CORR=TRUE, PROBT1=TRUE, DIC=TRUE, PROBD=TRUE,  probs= c(0.025, 0.50, 0.975), ...){        
# Note: I made some modifications. More modifications are needed to let start and end work well 
# for other values other than the default ones. 
#start <- burnin + 1
burnin <- start - 1
if(missing(end)) end   <- object$iter

if(start>=end) stop("The argument burnin should be smaller than the total number of iterations!")

lenProb    <- length(probs)
N          <- object$N
nTest      <- object$nTest
test.names <- object$test.names
nClass     <- object$nClass
nA         <- object$nA
nB         <- object$nB
patterns   <- object$patterns
bNzero     <- (nB!=0)

if(missing(D)){
DVec <- 1:nClass
names(DVec) <- paste("D",seq(1,nClass),   sep="")
D <- as.list(DVec)
}

# summarize prevalence
summary.prev <- t(rbind(apply(object$prev[start:end, ], 2, quantile,  probs= probs,na.rm = TRUE), 
      mean=apply(object$prev[start:end, ], 2, mean,na.rm = TRUE)))
rownames(summary.prev) <- paste("Class", as.character(1:nClass),  sep="")
# summarize parameter a
summary.a <- t(rbind(apply(object$a[start:end, ], 2, quantile,  probs= probs,na.rm = TRUE), 
      mean=apply(object$a[start:end, ], 2, mean,na.rm = TRUE)))
rownames(summary.a) <- paste("a", as.character(1:nA), sep="")
# summarize parameter b
if(nB>1){      
      summary.b <- t(rbind(apply(object$b[start:end, ], 2, quantile,  probs= probs,na.rm = TRUE), 
      mean=apply(object$b[start:end, ], 2, mean,na.rm = TRUE)))
      rownames(summary.b) <- paste("b", as.character(1:nB), sep="")
}else{
      if(nB==1){      
      summary.b <- t(rbind(quantile(object$b[start:end],   probs= probs,na.rm = TRUE), 
      mean=mean(object$b[start:end],na.rm = TRUE)))
#      rownames(b) <- paste("b", as.character(1:nB), sep="")                              
      }else{
      summary.b <- NULL      
      }
 }

    nPair <- nTest*(nTest-1)/2
    name.pair <-  rep(0, nPair)
    k <- 0
    for(i1 in (1:(nTest - 1))){
        for(i2 in ((i1+1):nTest)){       
         k <- k + 1
            name.pair[k] <- paste(test.names[i1], test.names[i2],  sep="&")              
        }
    }   # end of for(i1 in (1:(nTest - 1)))


AGREE=TRUE;
if(AGREE){
    prEgO <- t(apply(object$agree$pAgree[,, start:end]>object$agree$agree[,, start:end], c(1, 2), mean,na.rm = TRUE))
    colnames(prEgO) <- name.pair
    rownames(prEgO) <- paste("Class", 1:nClass,  sep="")

    str1 <- "list("
    for(i in 1:(nClass-1)) str1 <- paste(str1, "Class", i, "=NULL,", sep="")
    str1 <- paste(str1, "Class", nClass, "=NULL)", sep="")
    agree <- eval(parse(text=str1))                 
    pAgree <- agree 
    
    quantArray  <- apply(object$agree$pAgree[, , start:end], c(1, 2), quantile,  probs= probs, na.rm = TRUE)
    meanMat     <- apply(object$agree$pAgree[, , start:end], c(1, 2), mean,na.rm = TRUE)
    quantArray2 <- apply(object$agree$agree[, , start:end], c(1, 2), quantile,  probs= probs,na.rm = TRUE)
    meanMat2    <- apply(object$agree$agree[, , start:end], c(1, 2), mean,na.rm = TRUE)

    for(i in 1:nClass){
    pAgree[i] <- list(t(rbind(quantArray[,,i], mean=meanMat[,i])))
    agree[i]  <- list(t(rbind(quantArray2[,,i], mean=meanMat2[,i])))
    }
   agreelist=list(agree=agree, pAgree=pAgree);
   }else{
        prEgO = NULL;
        agreelist = NULL;
}   # end of if(AGREE);
    

sumy.gs <- list(call=object$call, nA=object$nA, nB=object$nB, nClass=object$nClass, bNzero=object$bNzero, nClass=object$nClass, init=object$init, prior.sdA=object$prior.sdA, sdA = object$sdA, prior.mA=object$prior.mA, mA=object$mA, 
    priorPrev=object$priorPrev, priorA=object$priorA, priorB=object$priorB, 	iter=object$iter, allCell=object$allCell, burnin=burnin, 
	summary.prev=summary.prev, summary.a=summary.a, summary.b=summary.b, 
    summary.prEgO = prEgO, summary.agree=agreelist)


PRED=TRUE;
if(PRED){
      pred <- t(rbind(apply(object$predmat[start:end, ], 2, quantile,  probs= probs, na.rm = TRUE),
                  mean = apply(object$predmat[start:end, ], 2, mean, na.rm = TRUE)))
      rownames(pred) <- patterns
      sumy.gs$summary.pred <- pred
    }


if(CORR){
    ptm = proc.time()
    cat("Calculating the correlation residuals...\n")
    try(sumy.gs$corr <- corr(object, start, end))
    print(proc.time() -ptm)
	pCorr <- t(rbind(apply(sumy.gs$corr$pCorr, 1, quantile,  probs= probs, na.rm = TRUE),
                  mean = apply(sumy.gs$corr$pCorr, 1, mean, na.rm = TRUE)))

           pCorr2 <- t(rbind(apply(-sumy.gs$corr$pCorr, 1, quantile,  probs= probs, na.rm = TRUE),
                  mean = apply(-sumy.gs$corr$pCorr, 1, mean, na.rm = TRUE)))

           corr <- matrix(1, nrow=ncol(pCorr), ncol=1)%*%(sumy.gs$corr$corr) 
           corr <- t(corr)
           rownames(pCorr) <-  name.pair
           resid  <- corr + pCorr2     
           rownames(resid) <-  name.pair 
           sumy.gs$summary.corr <- list(corr=sumy.gs$corr$corr, pCorr=pCorr, resid=resid)       
}

if(PROBT1){
    ptm = proc.time()
    cat("Calculating Pr(positive result|latent class)...\n")
    try(sumy.gs$probT1 <- probT1(object,  start,end))
	object$probT1 <- sumy.gs$probT1
    print(proc.time() -ptm)
        str1 <- "list("
        for(i in 1:(nClass-1)) str1 <- paste(str1, "Class", i, "=NULL,", sep="")
        str1 <- paste(str1, "Class", nClass, "=NULL)", sep="")
        probT1 <- eval(parse(text=str1))                    
            
        quantArray <- apply(sumy.gs$probT1, c(2, 1), quantile,  probs= probs, na.rm = TRUE)
        meanMat <- apply(sumy.gs$probT1, c(2, 1), mean, na.rm = TRUE)
    
        for(i in 1:nClass){ 
            mat <- t(rbind(quantArray[,,i], mean=meanMat[,i]))
            str.first <-  paste("P(", test.names, "+", sep="")
            str.last <- paste("|Class", i, ")", sep="")
            rownames(mat) <- paste(str.first, str.last, sep="")                         
            probT1[i] <- list(mat)
        }                   
        sumy.gs$summary.probT1 <- probT1
	
    }

#DIC=TRUE;
if(DIC){
   ptm = proc.time()
   cat("Calculating log-likelihood and DIC...\n")   
   try(redic  <- logLik(object,  start, end))
   sumy.gs$logLik  <- redic$logLik
   sumy.gs$Dbar    <- redic$dbar                
   sumy.gs$pD      <- redic$pd
   sumy.gs$DIC     <- redic$dic
   sumy.gs$sumy.logLik <- redic$sumy.loglik 
   print(proc.time() -ptm)
   }

if(PROBD){
    try(sumy.gs$summary.probD <- probD(object,  start,end, probs))
    }
    
  Dnames <- names(D)      
  ss.result <- list()
        for(i in 1:length(D)){              
            ss.out <- ss(object, start,end, D=as.numeric(unlist(D[i])), probs)                       
            str1 <- paste("ss.result$", Dnames[i], "<- list(test=as.numeric(unlist(D[i])), sensitivity=ss.out$sensitivity, specificity=ss.out$specificity)",  sep="")           
            eval(parse(text=str1))                                  
        }
  sumy.gs$summary.ss <- ss.result                     
  class(sumy.gs) <- c("summary.gibbs")
  return(sumy.gs)
}



#=================print.summary.gibbs=============================
print.summary.gibbs <- function(x,  ...){    
	cat("Lower limit and the upper limit of the prior distribution of a's\n", x$priorA$l, ",", x$priorA$u, "\n")
	cat("The priors for mA are\n", x$prior.mA, "\n")
	cat("The priors for sdA are \n", x$prior.sdA, "\n")
	cat("Lower limit and the upper limit of the prior distribution of b's \n", x$priorB$l, ",", x$priorB$u, "\n")
	cat("The priors for prevalence are \n", x$priorPrev, "\n")
	
    print(x$call)
    cat("Iteration: ", x$burnin+1, "--", x$iter, "\n")
    cat("Prevalence:\n")

    print(round(x$summary.prev, 3))
    cat("\n Parameter a's:\n")
    print(round(x$summary.a, 2))
    if(x$bNzero){
    cat("\n Parameter b's:\n")
    print(round(x$summary.b, 2))
    }

    if(!is.null(x$summary.corr)){
    cat("\n Observed Correlation:\n")
    print(round(x$summary.corr$corr, 3))
    cat("\n Expected Correlation:\n")
    print(round(x$summary.corr$pCorr, 3))
    cat("\n Correlation residuals:\n")
    print(round(x$summary.corr$resid, 3))
    }

    if(!is.null(x$summary.agree)){
    cat("\n Agreement in the expected:\n")
    for(i in 1:x$nClass){
    cat("Class", i, "\n")
    print(round(eval(parse(text=paste("x$summary.agree$pAgree$Class", i, sep=""))), 3))
    }



    cat("\n Agreement in the observed:\n")
    for(i in 1:x$nClass){
    cat("Class", i, "\n")
    print(round(eval(parse(text=paste("x$summary.agree$agree$Class", i, sep=""))), 3))
    }

    cat("\n P(E>O):\n")
    print(round(x$summary.prEgO, 3))
    }

    if(!is.null(x$summary.pred)){
    cat("\n Predicted values:\n")
    print(cbind(round(x$summary.pred, 0), Observed=x$allCell))
    }

    if(!is.null(x$summary.probT1)){
    cat("\n Pr(T+|Class):\n")
    for(i in 1:x$nClass){
    cat("Class", i, "\n")
    print(round(eval(parse(text=paste("x$summary.probT1$Class", i, sep=""))),3) )
    }
    }

    Dnames <- names(x$summary.ss)
    for(i in 1:length(x$summary.ss)){
        ss.out <- eval(parse(text=paste("x$summary.ss$", Dnames[i], sep="")))
        if(!is.null(ss.out$sensitivity)){
        cat("\n Sensitivities for", Dnames[i],"(class(es)", ss.out$test, "):\n")
        print(round(ss.out$sensitivity, 3))
        cat("\n Specificities for", Dnames[i],"(class(es)", ss.out$test, "):\n")
        print(round(ss.out$specificity, 3))
        }
    }


  #  if(!is.null(x$Dbar)){
  #  cat("\n Number of parameters: ", x$nA+x$nB+x$nClass-1,  "; DIC:",  round(x$DIC, 3),  "; \n logLik:",  round(x$sumy.logLik, 4), "\n Dbar:",  round(x$Dbar, 3),"; pD:",  round(x$pD, 3), "\n")
  #  }

   if(!is.null(x$logLik)){
    sumy.logLik <- quantile(x$logLik, c(0.025, 0.50, 0.975))   
    cat("\n Number of parameters: ", x$nA+x$nB+x$nClass-1, "; \n logLik:\n")
	print(sumy.logLik)
	cat("DIC:",  round(x$DIC, 3),   ";  Dbar:",  round(x$Dbar, 3),"; pD:",  round(x$pD, 3), "\n")
    }

    if(!is.null(x$summary.probD)){
    cat("\n Pr(Class+|Pattern):\n")
    for(i in 1:x$nClass){
    cat("Class", i, "\n")
    print(round(eval(parse(text=paste("x$summary.probD$Class", i, sep=""))), 2))
    }
    }
}







#-----------------------------------------------------------------------------
#gauQuad <- function(numGQ=20){
#   # default number of the gaussian quadrature values is 20        
#   loadSM <- library(statmod, logical.return=TRUE) 
#   if(loadSM){
#   out <- gauss.quad(numGQ,"hermite")
#   gq <- c(sqrt(2)*as.vector(out$nodes))                  # gaussian quadrature values 
#   w <- as.vector(out$weights)/sqrt(3.141593) 
#   }else{  
#   numGQ <- 20
#   gq <- c(-7.619049, -6.51059, -5.578739,-4.734581,-3.943967,-3.189015,-2.458664,-1.745247,-1.042945,-0.3469642,0.3469642,1.042945,1.745247,2.458664,3.189015,3.943967,4.734581,5.578739,6.51059,7.619049)
#   w  <- c(1.257801e-13,2.482062e-10,6.12749e-08,4.402121e-06,0.0001288263,0.001830103,0.01399784,0.06150637,0.1617393,0.2607930,0.2607930,0.1617393,0.06150637,0.01399784,0.001830103,0.0001288263,4.402121e-06,6.12749e-08,2.482062e-10,1.257801e-13)
#   }
#   list(numGQ=numGQ, gq=gq, w=w)
#}




#===================== plot.gibbs =================================
plot.gibbs <- function(x, option=1,ask=TRUE, ...){ 
if(class(x)!="gibbs") stop("x is an object of class gibbs!\n")  
# option=1: check the convergence; option=2: look at the agreement.

    oldask <- par("ask")
    par("ask"=ask)
    nClass <-  x$nClass
    J      <-  x$nTest
    nA     <-  x$nA
    nB     <-  x$nB 
    prev   <- rep(0, nClass)        
    bNzero <- (nB!=0)   
    niter  <- nrow(x$prev)
if(option==1){
for(i in 1:nClass){       
    plot(1:niter, x$prev[, i], xlab="Iteration Number", ylab=paste("Prevalence", i))    
    }
    for(i in 1:nA){
    plot(1:niter, x$a[, i], xlab="Iteration Number", ylab=paste("Parameter a", i, sep=""))  
    }
if(bNzero){
for(i in 1:nB){
    plot(1:niter, x$b[, i], xlab="Iteration Number", ylab=paste("Parameter b", i, sep=""))  
    }
    }
}

if(option==2){      
# nSubplot <- J*(J-1)/2
x.label <- "Predicted agreement"
y.label <- "Observed agreement"

for(n.i in 1:nClass){
    k <- 0
#   par(mfrow=parRC(nSubplot))
    for(i1 in (1:(J-1))){
    nSubplot <- J-i1;
    par(mfrow=parRC(nSubplot))
        for(i2 in ((i1+1):J)){       
            k <- k + 1          
            main.label <-  paste("L", n.i, ": T", i1, "&", i2, sep="")
            plot(x$agree$pAgree[k, n.i, ], x$agree$agree[k, n.i, ], ylim=c(min(x$agree$pAgree[k, n.i, ],x$agree$agree[k, n.i, ], na.rm=TRUE),max(x$agree$pAgree[k, n.i, ],x$agree$agree[k, n.i, ], na.rm=TRUE)), main=main.label, xlab=x.label, ylab=y.label)
            lines(x$agree$pAgree[k, n.i, ], x$agree$pAgree[k, n.i, ])
        }
    }
}   
} # end of if(option==2)
par("ask"=oldask)
}        
