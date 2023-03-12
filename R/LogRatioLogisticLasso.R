#' LogRatioLogisticLasso: Logistic log-ratio lasso regression
#'
#' @description Conduct logistic log-ratio lasso regression
#' @param x Covariate data matrix
#' @param y Binary outcome data vector
#' @param length.lambda Number of penalty parameters used in the path
#' @param lambda.min.ratio Ratio between the minimum and maximum choice of lambda. Default is `NULL`, where the ratio is chosen as 1e-2 if n < p and 1e-4 otherwise.
#' @param mu Value of penalty for the augmented Lagrangian
#' @param ncv Number of cross-validation runs. Use `NULL` if cross-validation is not wanted.
#' @param foldid A vector of fold indicator. Default is `NULL`.
#' @param step2 TRUE or FALSE, indicating whether a stepwise feature selection should be performed for features selected by the main lasso algorithm. Will only be performed if cross validation is enabled.
#' @param progress TRUE or FALSE, indicating whether printing progress bar as the algorithm runs.
#' @param plot TRUE or FALSE, indicating whether returning plots of model fitting.
#' @return A list with path-specific estimates (beta), path (lambda), and many others.
#' @author Teng Fei. Email: feit1@mskcc.org
#'
#' @import Rcpp RcppArmadillo ggplot2 reshape RcppProgress
#' @useDynLib LogRatioReg
#' @export

LogRatioLogisticLasso <- function(x,
                          y,
                          length.lambda=100,
                          lambda.min.ratio=NULL,
                          mu=1,
                          ncv=5,
                          foldid=NULL,
                          step2=FALSE,
                          progress=TRUE,
                          plot=TRUE){
  
  ptm <- proc.time()
  
  n <- length(y)
  p <- ncol(x)
  sfun = y-0.5
  lambda0 <- max(t(sfun) %*% x)/n
  
  if (is.null(lambda.min.ratio)) lambda.min.ratio = ifelse(n < p, 1e-01, 1e-02)
  
  lambda <- 10^(seq(log10(lambda0),log10(lambda0*lambda.min.ratio),length.out=length.lambda))
  
  if (progress) cat("Algorithm running for full dataset: \n")
  
  fullfit <- logistic_lasso_al(x,y,length.lambda,mu,100,lambda,progress)
  
  if (!is.null(colnames(x))){
    rownames(fullfit$beta) = colnames(x)
  }else{
    rownames(fullfit$beta) = 1:ncol(x)
  }
  
  if (!is.null(ncv)){
    
    cvmse <- matrix(NA,nrow=length.lambda,ncol=ncv)
    
    if (is.null(foldid)){
      labels <- caret::createFolds(factor(y),k=ncv,list=FALSE)
    }else{
      labels <- foldid
    }
    
    for (cv in 1:ncv){
      
      if (progress) cat(paste0("Algorithm running for cv dataset ",cv," out of ",ncv,": \n"))
      
      train.x <- x[labels!=cv,]
      train.y <- y[labels!=cv]
      test.x <- x[labels==cv,]
      test.y <- y[labels==cv]
      
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
                best.idx=best.idx,
                foldid=labels
    )
    
    if (plot){
      
      beta_nzero <- suppressWarnings(data.frame(reshape::melt(ret$beta[rowSums(ret$beta != 0) > 0,])))
      beta_nzero$lambda <- ret$lambda[beta_nzero$X2]
      beta_nzero$loglambda <- log(beta_nzero$lambda)
      
      lambda_count <- data.frame(loglambda = log(ret$lambda),
                                 count = colSums(ret$beta != 0))
      lambda_count <- lambda_count[seq(5,nrow(lambda_count),length.out=10),]
      
      top10feat <- sort(ret$beta[,length(ret$lambda)])[c(1:5,(p-4):p)]
      top10name <- names(top10feat)
      
      pcoef <- ggplot(beta_nzero, aes(x=loglambda,y=value,group=X1,color=as.factor(X1))) + 
        geom_line() + 
        scale_color_manual(values=rainbow(sum(rowSums(ret$beta != 0) > 0))) + 
        theme_bw() + 
        theme(legend.position = "none") +
        xlab("log(lambda)") +
        ylab("Coefficient") +
        geom_vline(xintercept=log(ret$lambda[ret$best.idx$idx.min]),linetype="dashed",color="darkgrey")+
        geom_vline(xintercept=log(ret$lambda[ret$best.idx$idx.1se]),linetype="dotted",color="darkgrey")+
        annotate("text",x=min(beta_nzero$loglambda)-2,y=top10feat,label=top10name,hjust=0)+
        annotate("text",x=lambda_count$loglambda,y=max(beta_nzero$value)+0.2,label=as.character(lambda_count$count))+
        ggtitle("Coefficients versus log(lambda)")
      
      mseplot <- data.frame(loglambda=log(ret$lambda),
                            mse=ret$cvmse.mean,
                            se=ret$cvmse.se,
                            mseaddse=ret$cvmse.mean+ret$cvmse.se,
                            mseminse=ret$cvmse.mean-ret$cvmse.se)
      
      pmse <- ggplot(mseplot, aes(x=loglambda, y=mse)) +
        geom_errorbar(aes(ymin=mseminse,ymax=mseaddse),color="grey")+
        geom_point(color="red")+
        theme_bw() +
        xlab("log(lambda)") +
        ylab("Mean-Squared Error")+
        geom_vline(xintercept=log(ret$lambda[ret$best.idx$idx.min]),linetype="dashed",color="darkgrey")+
        geom_vline(xintercept=log(ret$lambda[ret$best.idx$idx.1se]),linetype="dotted",color="darkgrey")+
        annotate("text",x=lambda_count$loglambda,y=max(mseplot$mseaddse)+0.05,label=as.character(lambda_count$count))+
        ggtitle("Cross-validated MSE versus log(lambda)")
      
      ret$pcoef <- pcoef
      ret$pmse <- pmse
      
    }
    
    if (step2){
      
      if (length(which(ret$best.beta$min.mse!=0))){
        idxs <- combn(which(ret$best.beta$min.mse!=0),2)
        x.select.min <- matrix(NA,nrow=n,ncol=ncol(idxs))
        for (k in 1:ncol(idxs)){
          x.select.min[,k] <- x[,idxs[1,k]] - x[,idxs[2,k]]
        }
        df_step2 <- data.frame(y=y,x=x.select.min)
        step2fit <- step(glm(y~.,data=df_step2,family=binomial),trace=0)
        vars <- as.numeric(sapply(names(step2fit$coefficients),function(x) strsplit(x,split = "[.]")[[1]][2]))
        
        selected <- idxs[,vars]
        for (k1 in 1:nrow(selected)){
          for (k2 in 1:ncol(selected)){
            selected[k1,k2] <- colnames(x)[as.numeric(selected[k1,k2])]
          }
        }
        ret$step2.feature.min = selected
        ret$step2fit.min <- step2fit
      }
      
      if (length(which(ret$best.beta$add.1se!=0))){
        idxs <- combn(which(ret$best.beta$add.1se!=0),2)
        x.select.min <- matrix(NA,nrow=n,ncol=ncol(idxs))
        for (k in 1:ncol(idxs)){
          x.select.min[,k] <- x[,idxs[1,k]] - x[,idxs[2,k]]
        }
        df_step2 <- data.frame(y=y,x=x.select.min)
        step2fit <- step(glm(y~.,data=df_step2,family=binomial),trace=0)
        vars <- as.numeric(sapply(names(step2fit$coefficients),function(x) strsplit(x,split = "[.]")[[1]][2]))
        
        selected <- idxs[,vars]
        for (k1 in 1:nrow(selected)){
          for (k2 in 1:ncol(selected)){
            selected[k1,k2] <- colnames(x)[as.numeric(selected[k1,k2])]
          }
        }
        ret$step2.feature.1se = selected
        ret$step2fit.1se <- step2fit
      }
      
    }
    
  }else{
    
    ret <- list(beta=fullfit$beta,
                beta0=fullfit$beta0,
                lambda=fullfit$lambda,
                loss=fullfit$loss,
                mse=fullfit$mse
    )
    
    if (plot){
      
      beta_nzero <- suppressWarnings(data.frame(reshape::melt(ret$beta[rowSums(ret$beta != 0) > 0,])))
      beta_nzero$lambda <- ret$lambda[beta_nzero$X2]
      beta_nzero$loglambda <- log(beta_nzero$lambda)
      
      lambda_count <- data.frame(loglambda = log(ret$lambda),
                                 count = colSums(ret$beta != 0))
      lambda_count <- lambda_count[seq(5,nrow(lambda_count),length.out=10),]
      
      top10feat <- sort(ret$beta[,length(ret$lambda)])[c(1:5,(p-4):p)]
      top10name <- names(top10feat)
      
      pcoef <- ggplot(beta_nzero, aes(x=loglambda,y=value,group=X1,color=as.factor(X1))) + 
        geom_line() + 
        scale_color_manual(values=rainbow(sum(rowSums(ret$beta != 0) > 0))) + 
        theme_bw() + 
        theme(legend.position = "none") +
        xlab("log(lambda)") +
        ylab("Coefficient") +
        annotate("text",x=min(beta_nzero$loglambda)-2,y=top10feat,label=top10name,hjust=0)+
        annotate("text",x=lambda_count$loglambda,y=max(beta_nzero$value)+0.2,label=as.character(lambda_count$count))+
        ggtitle("Coefficients versus log(lambda)")
      
      ret$pcoef <- pcoef
      
    }
    
  }
  
  ret$runtime <- proc.time() - ptm
  
  return(ret)
  
}