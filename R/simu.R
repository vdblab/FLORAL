#' Simulate data following log-ratio model
#'
#' @description Simulate a dataset from log-ratio model.
#' @param n An integer of sample size
#' @param p An integer of number of features (taxa).
#' @param model Type of models associated with outcome variable, can be "linear", "binomial", "cox", "finegray", or "timedep" (survival endpoint with time-dependent features).
#' @param weak Number of features with \code{weak} effect size.
#' @param strong Number of features with \code{strong} effect size.
#' @param weaksize Actual effect size for \code{weak} effect size. Must be positive.
#' @param strongsize Actual effect size for \code{strong} effect size. Must be positive.
#' @param pct.sparsity Percentage of zero counts for each sample.
#' @param rho Parameter controlling the correlated structure between taxa. Ranges between 0 and 1.
#' @param ncov Number of covariates that are not compositional features.
#' @param betacov Coefficients corresponding to the covariates that are not compositional features.
#' @param method Simulation schemes. Options are "manual" and "SparseDOSSA2", where the "manual" version is described in the preprint and "SparseDOSSA2" uses the default stool template from the package `SparseDOSSA2`. Note that only sample size and effect size are required for SparseDOSSA2 simulation.
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
#' @importFrom msm dpexp ppexp rpexp
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
                 timedep_slope=NULL,
                 timedep_cor=NULL,
                 ncov=0,
                 betacov=0,
                 intercept=FALSE,
                 method="manual"){
  
  true_set <- 1:(weak+strong)
  weak_idx <- 1:weak
  strong_idx <- (weak+1):(weak+strong)
  
  beta <- rep(NA,weak+strong)
  beta[weak_idx] <- rep(c(weaksize,-weaksize),weak/2)
  beta[strong_idx] <- rep(c(strongsize,-strongsize),strong/2)
  
  if (model %in% c("linear","binomial")){
    y <- rep(NA,length=n)
  }else if (model == "cox"){
    t <- rep(NA,length=n)
    d <- rep(NA,length=n)
  }else if (model == "finegray"){
    t <- t0 <- rep(NA,length=n)
    d <- rep(NA,length=n)
  }else if (model == "timedep"){
    m <- 10
    id.vect <- rep(1:n, each = m)
    n0 <- n
    n <- length(id.vect)
    
    if (is.null(timedep_cor)){
      timedep_cor <- 0.4
    }
    if (is.null(timedep_slope)){
      timedep_slope <- 0.5
    }else{
      timedep_slope <- timedep_slope
    }
  }
  
  x <- xobs <- matrix(NA,nrow=n,ncol=p)
  
  seqdep <- floor(runif(n,min=5000,max=50000))
  # highidx <- true_set
  
  ###############################################
  
  if (method == "manual"){
    
    if (model != "timedep"){
      
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
      
    }else if(model == "timedep"){
      
      sigma <- rho^(as.matrix(dist(1:(p-(weak+strong)))))
      diag(sigma) <- 1
      mu <- rep(0,p-(weak+strong))
      
      x0 <- mvtnorm::rmvnorm(n=n,mean=mu,sigma=sigma)
      x[,setdiff(1:p,true_set)] <- x0
      
      for (i in 1:n0){
        
        # Mu <- rep(c(rep(log(p),3),0,rep(log(p),2),0,log(p),rep(0,2)),m)
        # mu <- c(rep(log(p),3),0,rep(log(p),2),0,log(p),rep(0,2))
        # slopes <- outer(0:9,timedep_slope)
        # Mu <- t(mu + t(slopes))
        
        slopes <- rep(0:1,(weak+strong)/2)*timedep_slope
        Mu <- t(c(rep(log(p),3),0,rep(log(p),2),0,log(p),rep(0,2)) + t(outer(0:9,slopes)))
        
        sigma1 <- rho^(as.matrix(dist(1:(weak+strong))))
        Sigma <- (diag(m) %x% sigma1) + ((matrix(1,nrow=m,ncol=m) - diag(m)) %x% (diag(strong+weak)*timedep_cor))
        # Sigma <- Matrix::bdiag(replicate(m,sigma,simplify=FALSE))
        # sigma_offdiag <- diag(strong+weak)*0.8
        # mat_template <- matrix(1,nrow=m,ncol=m) - diag(m)
        
        x1 <- matrix(mvtnorm::rmvnorm(n=1,mean=Mu,sigma=Sigma),nrow=m,ncol=weak+strong,byrow=TRUE)
        x[id.vect==i,true_set] <- x1
        
        if (pct.sparsity > 0){
          # for (i in 1:n0){
          zeroidx <- sample(true_set,size=floor(length(true_set)*pct.sparsity))
          x[id.vect==i,zeroidx] <- -Inf
          # }
        }
        
      }
      
      if (pct.sparsity > 0){
        for (j in 1:n){
          zeroidx <- sample(setdiff(1:p,true_set),size=floor(p*pct.sparsity))
          x[j,zeroidx] <- -Inf
        }
      }
      
    }
    
    x = apply(x,2,function(y) exp(y)/rowSums(exp(x)))
    
    for (k in 1:n){
      xobs[k,] <- rmultinom(1,size=seqdep[k],prob=x[k,])
    }
    xcount = xobs
    colnames(xcount) <- paste0("taxa",1:p)
    
    for (k in 1:n){
      x[k,] <- rmultinom(1,size=1000000,prob=x[k,])
    }
    x = log(x+1)
    
  }else if (method == "SparseDOSSA2"){
    
    # sim <- SparseDOSSA2::SparseDOSSA2(template = "Stool",
    #                                   n_sample=n,
    #                                   median_read_depth = 25000,
    #                                   new_features=FALSE,
    #                                   verbose = FALSE)
    # xcount <- t(sim$simulated_data)
    # taxa <- colnames(xcount)
    # 
    # true_set <- which(taxa %in% c("k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Blautia|s__Ruminococcus_torques",
    #                               "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiaceae|g__Clostridium|s__Clostridium_leptum",
    #                               "k__Bacteria|p__Firmicutes|c__Negativicutes|o__Selenomonadales|f__Veillonellaceae|g__Veillonella|s__Veillonella_unclassified",
    #                               "k__Bacteria|p__Verrucomicrobia|c__Verrucomicrobiae|o__Verrucomicrobiales|f__Verrucomicrobiaceae|g__Akkermansia|s__Akkermansia_muciniphila",
    #                               "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_uniformis",
    #                               "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Porphyromonadaceae|g__Parabacteroides|s__Parabacteroides_merdae",
    #                               "k__Bacteria|p__Firmicutes|c__Negativicutes|o__Selenomonadales|f__Veillonellaceae|g__Veillonella|s__Veillonella_parvula",
    #                               "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Ruminococcus|s__Ruminococcus_bromii",
    #                               "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Roseburia|s__Roseburia_inulinivorans",
    #                               "k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Coriobacteriales|f__Coriobacteriaceae|g__Collinsella|s__Collinsella_aerofaciens"
    # ))
    # 
    # colnames(xcount) <- paste0("taxa",1:ncol(xcount))
    # rownames(xcount) <- NULL
    # x = t(sim$simulated_matrices$rel)
    # for (k in 1:nrow(xcount)){
    #   x[k,] <- rmultinom(1,size=1000000,prob=x[k,])
    # }
    # x = log(x+1)
    
  }
  
  betavec <- rep(0,ncol(xcount))
  betavec[true_set] <- beta
  xobs <- log(xcount+1)
  colnames(xobs) <- NULL
  
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
    
    ret <- list(xcount=xcount,x=xobs,y=y,y0=y0,beta=betavec,idx=true_set)
    
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
    
    ret <- list(xcount=xcount,x=xobs,y=y,beta=betavec,idx=true_set)
    
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
    
    ret <- list(xcount=xcount,x=xobs,t=t,d=d,beta=betavec,idx=true_set)
    
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
    
    ret <- list(xcount=xcount,x=xobs,t=t,d=d,beta=betavec,idx=true_set)
    
    if (ncov > 0) ret$xcov=xcov
    
  }else if (model=="timedep"){
    
    g.inv <- sqrt
    g <- function(x) {
      x^2
    }
    
    t.max <- 9
    t.min <- 0
    
    z.list <- list()
    for (i in 1:n0) {
      z <- x[id.vect == i,true_set]
      z.list[[i]] <- cbind(z, exp(z %*% beta))
    }
    
    k <- function(x, m, M, rates, t){
      ifelse(x <= m | x >= M, 0, dpexp(x, rates, t))
    }
    
    gen.y <- function(z,time_base,t.max,t.min){
      
      exp.etay <- z[,ncol(z)]
      
      t <- time_base
      t.diff <- (t[-1] - t[1:(length(t) - 1)])[-(length(t) - 1)]
      g.inv.t <- g.inv(t)
      g.inv.t.diff <- (g.inv(t[-1]) - g.inv(t[1:(length(t) - 1)]))[-(length(t) - 1)]
      
      g.inv.t.max <- g.inv(t.max)
      g.inv.t.min <- g.inv(t.min)
      
      x1 <- exp.etay
      d <- ppexp(g.inv.t.max, x1, g.inv.t) - ppexp(g.inv.t.min, x1,g.inv.t)
      M <- 1 / d
      r <- 60
      repeat{
        y <- rpexp(r, x1, g.inv.t)
        u <- runif(r)
        t <- M * ((k(y, g.inv.t.min, g.inv.t.max, x1, g.inv.t) / d /
                     dpexp(y, x1, g.inv.t)))
        y <- y[u <= t][1]
        if (!is.na(y)) break
      }
      y
    }
    
    y <- sapply(z.list, gen.y, time_base = 0:9, t.max=t.max, t.min=t.min)
    g.y <- g(y)
    
    # print(summary(g.y))
    
    # ct <- runif(n, min=5,max=6)
    # d <- ifelse(g.y < ct, 1, 0)
    # g.y <- ifelse(g.y < ct, g.y, ct)
    
    prop.cen <- 0.5
    d <- sample(0:1, n, replace = TRUE, prob = c(prop.cen, 1 - prop.cen))
    
    data <- NULL
    for (i in 1:n0) {
      id.temp <- rep(i, ceiling(g.y[i]))
      time.temp <- c(1:ceiling(g.y[i]))
      time0.temp <- 0:(ceiling(g.y[i]) - 1)
      d.temp <- c(rep(0, length(time.temp) - 1), d[i])
      
      if (ceiling(g.y[i]) > 9){
        z.temp <- xcount[id.vect==i,]
        
        if (length(id.temp) > nrow(z.temp)){
          
          z.temp <- rbind(z.temp,
                          # do.call("rbind", replicate(length(id.temp) - nrow(z.temp),
                          #                            z.temp[nrow(z.temp),],
                          #                            simplify = FALSE))
                          z.temp[nrow(z.temp),]
          )
          
          id.temp <- id.temp[c(1:10,length(id.temp))]
          time.temp <- time.temp[c(1:10,length(id.temp))]
          time0.temp <- time0.temp[c(1:11)]
          d.temp <- d.temp[c(1:10,length(id.temp))]
          
        }
        
      }else{
        z.temp <- xcount[id.vect==i,][1:(ceiling(g.y[i])),]
      }
      
      if (is.null(nrow(z.temp))){
        
        data.temp <- c(id.temp, time.temp, time0.temp, d.temp, z.temp)
        
      }else{
        
        data.temp <- cbind(id.temp, time.temp, time0.temp, d.temp, z.temp)
      }
      
      data <- rbind(data, data.temp)
    }
    
    colnames(data) <- c('id', 't', 't0', 'd', paste0("z",1:p))
    data <- data.frame(data)
    data_unique <- data |> 
      group_by(id) |> 
      filter(row_number() == n()) |> 
      ungroup()
    
    xcount <- data[,-c(1:4)]
    colnames(xcount) <- rownames(xcount) <- NULL
    colnames(xcount) <- paste0("taxa",1:p)
    
    xcount_baseline <- data_unique[,-c(1:4)]
    colnames(xcount_baseline) <- rownames(xcount_baseline) <- NULL
    colnames(xcount_baseline) <- paste0("taxa",1:p)
    
    ret <- list(xcount=xcount,
                xcount_baseline=xcount_baseline,
                data=data,
                data_unique=data_unique,
                beta=betavec,
                idx=true_set)
    
  }
  
  return(ret)
  
}
