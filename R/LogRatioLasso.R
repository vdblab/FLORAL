LogRatioLasso <- function(x,
                          y,
                          length.lambda=100,
                          mu=1,
                          family="gaussian",
                          ncv=5,
                          intercept=FALSE){
  
  y0 <- y
  y <- y - mean(y)
  lambda0 <- max(t(y) %*% x)/nrow(x)
  n <- length(y)
  p <- ncol(x)
  
  lambda.min.ratio = ifelse(n < p, 1e-08, 1e-08)
  lambda <- 10^(seq(log10(lambda0),log10(lambda0*lambda.min.ratio),length.out=length.lambda))
  
  fullfit <- linear_lasso_al(x,y,length.lambda,mu,100,lambda,FALSE)
  
  if (intercept){
    beta0 <- mean(y0) - colMeans(x %*% fullfit$beta)
  }
  
  cvmse <- matrix(NA,nrow=length.lambda,ncol=ncv)
  
  if (!is.null(ncv)){
    
    labels <- caret::createFolds(1:nrow(x),k=ncv)
    
    for (cv in 1:ncv){
      
      train.x <- x[-labels[[cv]],]
      train.y <- y[-labels[[cv]]]
      test.x <- x[labels[[cv]],]
      test.y <- y[labels[[cv]]]
      
      cvfit <- linear_lasso_al(train.x,train.y,length.lambda,mu,100,lambda,FALSE)
      
      cvmse[,cv] <- apply(test.x %*% cvfit$beta,2,function(x) sum((test.y-x)^2)/length(test.y))
      
    }
    
  }
  
  mean.cvmse <- rowMeans(cvmse)
  se.cvmse <- apply(cvmse,1,sd)
  
  idx.min <- which.min(mean.cvmse)
  se.min <- se.cvmse[idx.min]
  idx.1se <- max(which(mean.cvmse > mean.cvmse[idx.min] + se.min & 1:length.lambda < idx.min))
  
  best.beta <- list(min.mse = fullfit$beta[,idx.min],
                    add.1se = fullfit$beta[,idx.1se])
  
  best.idx <- list(idx.min = idx.min,
                   idx.1se = idx.1se)
  
  ret <- list(beta=fullfit$beta,
              lambda=fullfit$lambda,
              loss=fullfit$loss,
              mse=fullfit$mse,
              cvmse.mean=mean.cvmse,
              cvmse.se=se.cvmse,
              best.beta=best.beta,
              best.idx=best.idx
              )
  
  if(intercept) ret$intercept = beta0
  
  return(ret)
  
}