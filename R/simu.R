#' Simulate data following log-ratio model
#'
#' @description Simulate a dataset from log-ratio model.
#' @param n An integer of sample size
#' @param p An integer of number of features (taxa).
#' @param model Type of models associated with outcome variable, can be "linear", "binomial", "cox", or "finegray".
#' @param weak Number of features with \code{weak} effect size.
#' @param strong Number of features with \code{strong} effect size.
#' @param weaksize Actual effect size for \code{weak} effect size. Must be positive.
#' @param strongsize Actual effect size for \code{strong} effect size. Must be positive.
#' @param pct.sparsity Percentage of zero counts for each sample.
#' @param rho Parameter controlling the correlated structure between taxa. Ranges between 0 and 1.
#' @param ncov Number of covariates that are not compositional features.
#' @param betacov Coefficients corresponding to the covariates that are not compositional features.
#' @param intercept Boolean. If TRUE, then a random intercept will be generated in the model. Only works for \code{linear} or \code{binomial} models.
#' @return A list with simulated count matrix \code{xcount}, log1p-transformed count matrix \code{x}, outcome (continuous \code{y}, continuous centered \code{y0}, binary \code{y}, or survival \code{t}, \code{d}), true coefficient vector \code{beta}, list of non-zero features \code{idx}, value of intercept \code{intercept} (if applicable).
#' @author Teng Fei. Email: feit1@mskcc.org
#'
#' @examples 
#' 
#' set.seed(23420)
#' dat <- simu(n=50,p=30,model="linear")
#' 
#' @import ggplot2 survival glmnet dplyr
#' @importFrom survcomp concordance.index
#' @importFrom reshape melt
#' @importFrom utils combn
#' @importFrom grDevices rainbow
#' @importFrom caret createFolds
#' @importFrom stats dist rbinom rexp rmultinom rnorm runif sd step glm binomial gaussian na.omit
#' @useDynLib FLORAL
#' @export

simu <- function(n = 100, 
                 p = 200, 
                 model = "linear",
                 weak = 4,
                 strong = 6,
                 weaksize = 0.125,
                 strongsize = 0.25,
                 pct.sparsity = 0.5,
                 rho=0,
                 ncov=0,
                 betacov=0,
                 intercept=FALSE){
  
  true_set <- 1:(weak+strong)
  weak_idx <- 1:weak
  strong_idx <- (weak+1):(weak+strong)
  
  beta <- rep(NA,weak+strong)
  beta[weak_idx] <- rep(c(weaksize,-weaksize),weak/2)
  beta[strong_idx] <- rep(c(strongsize,-strongsize),strong/2)
  
  x <- xobs <- matrix(NA,nrow=n,ncol=p)
  
  if (model %in% c("linear","binomial")){
    y <- rep(NA,length=n)
  }else if (model == "cox"){
    t <- rep(NA,length=n)
    d <- rep(NA,length=n)
  }else if (model == "finegray"){
    t <- t0 <- rep(NA,length=n)
    d <- rep(NA,length=n)
  }
  
  seqdep <- floor(runif(n,min=5000,max=50000))
  # highidx <- true_set
  
  sigma <- rho^(as.matrix(dist(1:p)))
  diag(sigma) <- c(rep(log(p)/2,3),1,rep(log(p)/2,2),1,log(p)/2,rep(1,p-8))
  mu <- c(rep(log(p),3),0,rep(log(p),2),0,log(p),rep(0,p-8))
  
  x <- mvtnorm::rmvnorm(n=n,mean=mu,sigma=sigma)
  
  if (pct.sparsity > 0){
    for (i in 1:n){
      zeroidx <- sample(1:p,size=floor(p*pct.sparsity))
      x[i,zeroidx] <- -Inf
    }
  }
  x = apply(x,2,function(y) exp(y)/rowSums(exp(x)))
  
  for (k in 1:n){
    xobs[k,] <- rmultinom(1,size=seqdep[k],prob=x[k,])
  }
  xobs = log(xobs+1)
  
  for (k in 1:n){
    x[k,] <- rmultinom(1,size=1000000,prob=x[k,])
  }
  xcount = x
  colnames(xcount) <- paste0("taxa",1:p)
  x = log(x+1)
  
  if (ncov > 0){
    xcov <- mvtnorm::rmvnorm(n=n,mean=rep(0,ncov))
    colnames(xcov) <- paste0("cov",1:ncov)
  }
  
  if (model == "linear"){
    
    y <- x[,true_set] %*% beta + rnorm(n,mean=0,sd=1)
    
    if(ncov > 0) y <- y + xcov %*% betacov
    
    if(intercept) {
      
      intcpt <- rnorm(1,mean=1,sd=1)
      y <- y + intcpt
      
    }
    y0 = y - mean(y)
    
    ret <- list(xcount=xcount,x=xobs,y=y,y0=y0,beta=c(beta,rep(0,p-weak-strong)),idx=true_set)
    
    if (intercept) ret$intercept=intcpt
    
    if (ncov > 0) ret$xcov=xcov
    
  }else if(model == "binomial"){
    
    eta <- x[,true_set] %*% beta
    
    if (ncov > 0) eta <- eta + xcov %*% betacov
    
    if(intercept) {
      
      intcpt <- rnorm(1,mean=1,sd=1)
      eta <- eta + intcpt
      
    }
    
    prob <- exp(eta)/(1+exp(eta))
    
    for (i in 1:n){
      y[i] <- rbinom(1,1,prob=prob[i])
    }
    
    ret <- list(xcount=xcount,x=xobs,y=y,beta=c(beta,rep(0,p-weak-strong)),idx=true_set)
    
    if (intercept) ret$intercept=intcpt
    
    if (ncov > 0) ret$xcov=xcov
    
  }else if(model == "cox"){
    
    eta <- x[,true_set] %*% beta
    
    if (ncov > 0) eta <- eta + xcov %*% betacov
    
    lambda <- exp(eta)
    
    for(i in 1:n){
      U <- runif(1,min=0,max=1)
      t0 <- log(1-(log(1-U))/(0.1*lambda[i]))
      c0 <- min(rexp(1,rate=0.1),runif(1,min=5,max=6))
      t[i] <- min(t0,c0)
      d[i] <- as.numeric(I(t0 <= c0))
    }
    
    ret <- list(xcount=xcount,x=xobs,t=t,d=d,beta=c(beta,rep(0,p-weak-strong)),idx=true_set)
    
    if (ncov > 0) ret$xcov=xcov
    
  }else if(model == "finegray"){
    
    eta <- x[,true_set] %*% beta
    
    if (ncov > 0) eta <- eta + xcov %*% betacov
    
    p.cif = 0.66
    lambda <- exp(eta)
    cl=0.19
    cu=10
    
    P1 <- 1-(1-p.cif)^(lambda)
    epsilon <- rep(0,n)
    for (i in 1:n){
      epsilon[i] <- 2 - rbinom(1,1,P1[i])
    }
    
    #generate the event time based on the type of outcome
    t0 <- rep(0,n)
    u <- runif(n)
    for (i in 1:n){
      if (epsilon[i] == 1){
        t0[i] <- -log(1 - (1 - (1-u[i]*(1-(1-p.cif)^lambda[i]))^(1/(lambda[i]+0.001)))/p.cif)
      }
      if (epsilon[i] == 2){
        t0[i] <- -log((1-u[i])^(1/(lambda[i]+0.001)))
      }
    }
    
    #generate censoring time
    c <- runif(n,cl,cu)
    #observed time
    t <- ifelse(t0 == Inf, c ,t0*I(t0<=c) + c*I(t0>c))
    # outcome
    d <- ifelse(t0 == Inf, 0, 0*I(t == c) + epsilon*I(t < c))
    
    ret <- list(xcount=xcount,x=xobs,t=t,d=d,beta=c(beta,rep(0,p-weak-strong)),idx=true_set)
    
    if (ncov > 0) ret$xcov=xcov
    
  }
  
  return(ret)
}
