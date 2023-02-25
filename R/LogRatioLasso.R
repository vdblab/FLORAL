#' LogRatioLasso: Linear log-ratio lasso regression
#'
#' @description Conduct linear log-ratio lasso regression
#' @param x Covariate data matrix
#' @param y Continuous outcome data vector
#' @param length.lambda Number of penalty parameters used in the path
#' @param mu Value of penalty for the augmented Lagrangian
#' @param ncv Number of cross-validation runs. Use `NULL` if cross-validation is not wanted.
#' @param intercept TRUE or FALSE, indicating whether an intercept should be estimated.
#' @param progress TRUE or FALSE, indicating whether printing progress bar as the algorithm runs.
#' @return A list with path-specific estimates (beta), path (lambda), and many others.
#' @author Teng Fei. Email: feit1@mskcc.org
#' 
#' @import caret Rcpp RcppArmadillo
#' @useDynLib LogRatioReg
#' @export
LogRatioLasso <- function(x,
                          y,
                          length.lambda=100,
                          mu=1,
                          ncv=5,
                          intercept=FALSE,
                          progress=TRUE){
  
  ptm <- proc.time()
  
  y0 <- y
  y <- y - mean(y)
  n <- length(y)
  p <- ncol(x)
  lambda0 <- max(t(y) %*% x)/n
  
  lambda.min.ratio = ifelse(n < p, 1e-08, 1e-08)
  lambda <- 10^(seq(log10(lambda0),log10(lambda0*lambda.min.ratio),length.out=length.lambda))
  
  if (progress) cat("Algorithm running for full dataset: \n")
  
  fullfit <- linear_lasso_al(x,y,length.lambda,mu,100,lambda,FALSE,progress)
  
  if (!is.null(colnames(x))) rownames(fullfit$beta) = colnames(x)
  
  if (intercept){
    beta0 <- mean(y0) - colMeans(x %*% fullfit$beta)
  }
  
  if (!is.null(ncv)){
    
    cvmse <- matrix(NA,nrow=length.lambda,ncol=ncv)
    
    labels <- caret::createFolds(1:nrow(x),k=ncv)
    
    for (cv in 1:ncv){
      
      if (progress) cat(paste0("Algorithm running for cv dataset ",cv," out of ",ncv,": \n"))
      
      train.x <- x[-labels[[cv]],]
      train.y <- y[-labels[[cv]]]
      test.x <- x[labels[[cv]],]
      test.y <- y[labels[[cv]]]
      
      cvfit <- linear_lasso_al(train.x,train.y,length.lambda,mu,100,lambda,FALSE,progress)
      
      cvmse[,cv] <- apply(test.x %*% cvfit$beta,2,function(x) sum((test.y-x)^2)/length(test.y))
      
    }
    
    
    mean.cvmse <- rowMeans(cvmse)
    se.cvmse <- apply(cvmse,1,sd)
    
    idx.min <- which.min(mean.cvmse)
    se.min <- se.cvmse[idx.min]
    idx.1se <- suppressWarnings(max(which(mean.cvmse > mean.cvmse[idx.min] + se.min & 1:length.lambda < idx.min)))
    if (idx.1se == -Inf) idx.1se = 1
    
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
    
  }else{
    
    ret <- list(beta=fullfit$beta,
                lambda=fullfit$lambda,
                loss=fullfit$loss,
                mse=fullfit$mse
    )
    
  }
  
  if(intercept) ret$intercept = beta0
  
  ret$runtime = proc.time() - ptm
  
  return(ret)
  
}