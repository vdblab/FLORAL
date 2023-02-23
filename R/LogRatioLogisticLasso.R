LogRatioLogisticLasso <- function(x,
                          y,
                          length.lambda=100,
                          mu=1,
                          ncv=5){
  
  ptm <- proc.time()
  
  n <- length(y)
  p <- ncol(x)
  sfun = y-0.5
  lambda0 <- max(t(sfun) %*% x)/n
  
  lambda.min.ratio = ifelse(n < p, 1e-02, 1e-04)
  lambda <- 10^(seq(log10(lambda0),log10(lambda0*lambda.min.ratio),length.out=length.lambda))
  
  fullfit <- logistic_lasso_al(x,y,length.lambda,mu,100,lambda)
  
  cvmse <- matrix(NA,nrow=length.lambda,ncol=ncv)
  
  if (!is.null(ncv)){
    
    labels <- caret::createFolds(factor(y),k=ncv)
    
    for (cv in 1:ncv){
      
      train.x <- x[-labels[[cv]],]
      train.y <- y[-labels[[cv]]]
      test.x <- x[labels[[cv]],]
      test.y <- y[labels[[cv]]]
      
      cvfit <- logistic_lasso_al(train.x,train.y,length.lambda,mu,100,lambda)
      
      cvmse[,cv] <- apply(cbind(1,test.x) %*% rbind(t(cvfit$beta0),cvfit$beta),2,function(x) sum((test.y- exp(x)/(1+exp(x)))^2)/length(test.y))
      
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