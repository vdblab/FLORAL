#' LogRatioCoxLasso: Cox log-ratio lasso regression
#'
#' @description Conduct Cox proportional hazards log-ratio lasso regression
#' @param x Covariate data matrix
#' @param y A `survival::Surv` object
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
  
  fullfit <- cox_lasso_al(x,t,d,tj,length.lambda,mu,100,lambda,devnull)
  lidx <- which(fullfit$loglik != 0 | !is.nan(fullfit$loglik))
  
  dev <- -2*fullfit$loglik[lidx]
  
  if (!is.null(colnames(x))) rownames(fullfit$beta) = colnames(x)
  
  cvdev <- matrix(0,nrow=length(lidx),ncol=ncv)
  cvdevnull <- rep(0,ncol=ncv)
  
  if (!is.null(ncv)){
    
    # labels <- caret::createFolds(factor(d),k=ncv)
    labels <- coxsplity(as.matrix(y),ncv)
    
    for (cv in 1:ncv){
      
      # train.x <- x[-labels[[cv]],]
      # train.d <- d[-labels[[cv]]]
      # train.t <- t[-labels[[cv]]]
      # test.x <- x[labels[[cv]],]
      # test.d <- d[labels[[cv]]]
      # test.t <- t[labels[[cv]]]
      
      train.x <- x[labels!=cv,]
      train.d <- d[labels!=cv]
      train.t <- t[labels!=cv]
      test.x <- x[labels==cv,]
      test.d <- d[labels==cv]
      test.t <- t[labels==cv]
      
      cv.devnull <- 0
      train.tj <- sort(train.t[train.d==1])
      for (j in 1:length(train.tj)){
        cv.devnull <- cv.devnull + log(sum(train.t >= train.tj[j]))
      }
      cv.devnull <- 2*cv.devnull
      
      cvfit <- cox_lasso_al(train.x,train.t,train.d,train.tj,length(lidx),mu,100,lambda[lidx],cv.devnull)
      
      cv.devnull <- 0
      loglik <- rep(0,length(lidx))
      linprod <- test.x %*% cvfit$beta
      
      test.tj <- sort(test.t[test.d==1])
      for (j in 1:length(test.tj)){
        cv.devnull <- cv.devnull + log(sum(test.t >= test.tj[j]))
        if (sum(test.t >= test.tj[j]) > 1){
          cvdev[,cv] <- cvdev[,cv] + linprod[test.t == test.tj[j],] - 
            log(colSums(exp(linprod[test.t >= test.tj[j],])) + 1e-8)
        }else if (sum(test.t >= test.tj[j]) == 1){
          cvdev[,cv] <- cvdev[,cv] + linprod[test.t == test.tj[j],] - 
            log(sum(exp(linprod[test.t >= test.tj[j],])) + 1e-8)
        }
        
        
        #- log(accu(link(widx))+1e-8)
      }
      cvdevnull[cv] <- 2*cv.devnull
      cvdev[,cv] <- -2*cvdev[,cv]
      # cvmse[,cv] <- apply(cbind(1,test.x) %*% rbind(t(cvfit$beta0),cvfit$beta),2,function(x) sum((test.y - exp(x)/(1+exp(x)))^2)/length(test.y))
      
    }
    
  }
  
  mean.cvdev <- rowMeans(cvdev)
  se.cvdev <- apply(cvdev,1,sd)
  
  idx.min <- which.min(mean.cvdev)
  se.min <- se.cvdev[idx.min]
  idx.1se <- max(which(mean.cvdev > mean.cvdev[idx.min] + se.min & 1:length(lidx) < idx.min))
  if (idx.1se == -Inf) idx.1se = 1
  
  best.beta <- list(min.mse = fullfit$beta[,idx.min],
                    add.1se = fullfit$beta[,idx.1se])
  
  best.idx <- list(idx.min = idx.min,
                   idx.1se = idx.1se)
  
  ret <- list(beta=fullfit$beta[,lidx],
              lambda=fullfit$lambda[lidx],
              loss=fullfit$loss[lidx],
              mse=fullfit$mse[lidx],
              cvdev.mean=mean.cvdev,
              cvdev.se=se.cvdev,
              best.beta=best.beta,
              best.idx=best.idx
  )
  
  ret$runtime <- proc.time() - ptm
  
  return(ret)
  
}