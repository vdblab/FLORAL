#' LogRatioCoxLasso: Cox log-ratio lasso regression
#'
#' @description Conduct Cox proportional hazards log-ratio lasso regression
#' @param x Covariate data matrix
#' @param y A `Surv` object
#' @param length.lambda Number of penalty parameters used in the path
#' @param mu Value of penalty for the augmented Lagrangian
#' @param ncv Number of cross-validation runs. Use `NULL` if cross-validation is not wanted.
#' @return A list with path-specific estimates (beta), path (lambda), and many others.
#' @author Teng Fei. Email: feit1@mskcc.org
#'
#' @import survival caret Rcpp RcppArmadillo
#' @useDynLib LogRatioReg
#' @export

LogRatioCoxLasso <- function(x,
                             y,
                             length.lambda=100,
                             mu=1,
                             ncv=5){
  
  ptm <- proc.time()
  
  t <- y[,1]
  d <- y[,2]
  n <- nrow(y)
  p <- ncol(x)
  
  tj <- sort(t[d==1])
  devnull <- log(sum(t >= tj[1]))
  denomj <- rep(NA,length(tj))
  denomj[1] <- sum(t >= tj[1])
  for (j in 2:length(tj)){
    denomj[j] <- denomj[j-1] - sum(t >= tj[j-1] & t < tj[j])
    devnull <- devnull + log(sum(t >= tj[j]))
  }
  devnull <- 2*devnull
  
  expect <- sapply(t,function(x) sum(1/denomj[which(tj <= x)]))
  sfun = d - expect
  # hfun <- expect - sapply(t,function(x) sum(1/denomj[which(tj <= x)]^2)) + 1e-8 # hessian
  # z <- sfun/hfun
  lambda0 <- max(t(sfun) %*% x)/n
  
  lambda.min.ratio = ifelse(n < p, 1e-02, 1e-04)
  lambda <- 10^(seq(log10(lambda0),log10(lambda0*lambda.min.ratio),length.out=length.lambda))
  
  fullfit <- cox_lasso_al(x,t,d,tj,length.lambda,mu,100,lambda)
  dev <- -2*fullfit$loglik
  
  if (!is.null(colnames(x))) rownames(fullfit$beta) = colnames(x)
  
  cvdevdiff <- matrix(0,nrow=length.lambda,ncol=ncv)
  cvdevnull <- rep(0,ncol=ncv)
  
  if (!is.null(ncv)){
    
    labels <- caret::createFolds(factor(d),k=ncv)
    
    for (cv in 1:ncv){
      
      train.x <- x[-labels[[cv]],]
      train.d <- d[-labels[[cv]]]
      train.t <- t[-labels[[cv]]]
      test.x <- x[labels[[cv]],]
      test.d <- d[labels[[cv]]]
      test.t <- t[labels[[cv]]]
      
      train.tj <- sort(train.t[train.d==1])
      
      cvfit <- cox_lasso_al(train.x,train.t,train.d,train.tj,length.lambda,mu,100,lambda)
      
      cv.devnull <- 0
      loglik <- rep(0,length.lambda)
      linprod <- test.x %*% cvfit$beta
      
      test.tj <- sort(test.t[test.d==1])
      for (j in 1:length(test.tj)){
        cv.devnull <- cv.devnull + log(sum(test.t >= test.tj[j]))
        cvdevdiff[,cv] <- cvdevdiff[,cv] + linprod[test.t == test.tj[j],] - 
          log(colSums(exp(linprod[test.t >= test.tj[j],])) + 1e-8)
        
        #- log(accu(link(widx))+1e-8)
      }
      cvdevnull[cv] <- 2*cv.devnull
      cvdevdiff[,cv] <- -2*cvdevdiff[,cv]
      # cvmse[,cv] <- apply(cbind(1,test.x) %*% rbind(t(cvfit$beta0),cvfit$beta),2,function(x) sum((test.y - exp(x)/(1+exp(x)))^2)/length(test.y))
      
    }
    
  }
  
  mean.cvmse <- rowMeans(cvmse)
  se.cvmse <- apply(cvmse,1,sd)
  
  idx.min <- which.min(mean.cvmse)
  se.min <- se.cvmse[idx.min]
  idx.1se <- max(which(mean.cvmse > mean.cvmse[idx.min] + se.min & 1:length.lambda < idx.min))
  if (idx.1se == -Inf) idx.1se = 1
  
  best.beta0 <- list(min.mse = fullfit$beta0[idx.min],
                     add.1se = fullfit$beta0[idx.1se])
  
  best.beta <- list(min.mse = fullfit$beta[,idx.min],
                    add.1se = fullfit$beta[,idx.1se])
  
  best.idx <- list(idx.min = idx.min,
                   idx.1se = idx.1se)
  
  ret <- list(beta=fullfit$beta,
              beta0=fullfit$beta0,
              lambda=fullfit$lambda,
              loss=fullfit$loss,
              mse=fullfit$mse,
              cvmse.mean=mean.cvmse,
              cvmse.se=se.cvmse,
              best.beta=best.beta,
              best.beta0=best.beta0,
              best.idx=best.idx
  )
  
  ret$runtime <- proc.time() - ptm
  
  return(ret)
  
}