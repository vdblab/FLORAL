#' LogRatioLogisticLasso: Logistic log-ratio lasso regression
#'
#' @description Conduct logistic log-ratio lasso regression
#' @param x Covariate data matrix
#' @param y Binary outcome data vector
#' @param length.lambda Number of penalty parameters used in the path
#' @param mu Value of penalty for the augmented Lagrangian
#' @param ncv Number of cross-validation runs. Use `NULL` if cross-validation is not wanted.
#' @param progress TRUE or FALSE, indicating whether printing progress bar as the algorithm runs.
#' @return A list with path-specific estimates (beta), path (lambda), and many others.
#' @author Teng Fei. Email: feit1@mskcc.org
#'
#' @import caret Rcpp RcppArmadillo
#' @useDynLib LogRatioReg
#' @export

LogRatioLogisticLasso <- function(x,
                          y,
                          length.lambda=100,
                          mu=1,
                          ncv=5,
                          progress=TRUE){
  
  ptm <- proc.time()
  
  n <- length(y)
  p <- ncol(x)
  sfun = y-0.5
  lambda0 <- max(t(sfun) %*% x)/n
  
  lambda.min.ratio = ifelse(n < p, 1e-02, 1e-04)
  lambda <- 10^(seq(log10(lambda0),log10(lambda0*lambda.min.ratio),length.out=length.lambda))
  
  if (progress) cat("Algorithm running for full dataset: \n")
  
  fullfit <- logistic_lasso_al(x,y,length.lambda,mu,100,lambda,progress)
  
  if (!is.null(colnames(x))) rownames(fullfit$beta) = colnames(x)
  
  if (!is.null(ncv)){
    
    cvmse <- matrix(NA,nrow=length.lambda,ncol=ncv)
    
    labels <- caret::createFolds(factor(y),k=ncv)
    
    for (cv in 1:ncv){
      
      if (progress) cat(paste0("Algorithm running for cv dataset ",cv," out of ",ncv,": \n"))
      
      train.x <- x[-labels[[cv]],]
      train.y <- y[-labels[[cv]]]
      test.x <- x[labels[[cv]],]
      test.y <- y[labels[[cv]]]
      
      cvfit <- logistic_lasso_al(train.x,train.y,length.lambda,mu,100,lambda,progress)
      
      cvmse[,cv] <- apply(cbind(1,test.x) %*% rbind(t(cvfit$beta0),cvfit$beta),2,function(x) sum((test.y- exp(x)/(1+exp(x)))^2)/length(test.y))
      
    }
    
    mean.cvmse <- rowMeans(cvmse)
    se.cvmse <- apply(cvmse,1,sd)
    
    idx.min <- which.min(mean.cvmse)
    se.min <- se.cvmse[idx.min]
    idx.1se <- suppressWarnings(max(which(mean.cvmse > mean.cvmse[idx.min] + se.min & 1:length.lambda < idx.min)))
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
    
  }else{
    
    ret <- list(beta=fullfit$beta,
                beta0=fullfit$beta0,
                lambda=fullfit$lambda,
                loss=fullfit$loss,
                mse=fullfit$mse
    )
    
  }
  
  ret$runtime <- proc.time() - ptm
  
  return(ret)
  
}