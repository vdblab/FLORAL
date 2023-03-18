#' LogRatioCoxLasso: Cox log-ratio lasso regression
#'
#' @description Conduct Cox proportional hazards log-ratio lasso regression
#' @param x Covariate data matrix
#' @param y A `survival::Surv` object
#' @param length.lambda Number of penalty parameters used in the path
#' @param lambda.min.ratio Ratio between the minimum and maximum choice of lambda. Default is `NULL`, where the ratio is chosen as 1e-2 if n < p and 1e-4 otherwise.
#' @param mu Value of penalty for the augmented Lagrangian
#' @param ncv Number of cross-validation runs. Use `NULL` if cross-validation is not wanted.
#' @param foldid A vector of fold indicator. Default is `NULL`.
#' @param step2 TRUE or FALSE, indicating whether a stepwise feature selection should be performed for features selected by the main lasso algorithm. Will only be performed if cross validation is enabled.
#' @param progress TRUE or FALSE, indicating whether printing progress bar as the algorithm runs.
#' @param plot TRUE or FALSE, indicating whether returning plots of model fitting.
#' @param mcv Metric for cross validation prediction assessment. Default is `Deviance`. An alternative option is `Cindex`.
#' @return A list with path-specific estimates (beta), path (lambda), and many others.
#' @author Teng Fei. Email: feit1@mskcc.org
#'
#' @import survival Rcpp RcppArmadillo ggplot2 reshape RcppProgress
#' @importFrom survcomp concordance.index
#' @useDynLib LogRatioReg
#' @export

LogRatioCoxLasso <- function(x,
                             y,
                             length.lambda=100,
                             lambda.min.ratio=NULL,
                             mu=1,
                             ncv=5,
                             foldid=NULL,
                             step2=FALSE,
                             progress=TRUE,
                             plot=TRUE,
                             mcv="Deviance"){
  
  ptm <- proc.time()
  
  t <- y[,1]
  d <- y[,2]
  n <- nrow(y)
  p <- ncol(x)
  if (length(unique(t)) < n)  t <- t+runif(n,min=0,max=1e-4) # break ties
  
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
  lambda0 <- 2*max(t(sfun) %*% x)/n
  
  if (is.null(lambda.min.ratio)) lambda.min.ratio = ifelse(n < p, 1e-02, 1e-02)
  lambda <- 10^(seq(log10(lambda0),log10(lambda0*lambda.min.ratio),length.out=length.lambda))
  
  if (progress) cat("Algorithm running for full dataset: \n")
  
  fullfit <- cox_lasso_al(x,t,d,tj,length.lambda,mu,100,lambda,devnull,progress)
  lidx <- which(fullfit$loglik != 0 | !is.nan(fullfit$loglik))
  
  dev <- -2*fullfit$loglik[lidx]
  
  if (!is.null(colnames(x))){
    rownames(fullfit$beta) = colnames(x)
  }else{
    rownames(fullfit$beta) = 1:ncol(x)
  }
  
  if (!is.null(ncv)){
    
    cvdev <- matrix(0,nrow=length(lidx),ncol=ncv)
    cvdevnull <- rep(0,ncol=ncv)
    
    if (is.null(foldid)){
      labels <- coxsplity(as.matrix(y),ncv)
    }else{
      labels <- foldid
    }
    
    # labels <- caret::createFolds(factor(d),k=ncv)
    
    for (cv in 1:ncv){
      
      if (progress) cat(paste0("Algorithm running for cv dataset ",cv," out of ",ncv,": \n"))
      
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
      
      cvfit <- cox_lasso_al(train.x,train.t,train.d,train.tj,length(lidx),mu,100,lambda[lidx],cv.devnull,progress)
      
      cv.devnull <- 0
      loglik <- rep(0,length(lidx))
      linprod <- test.x %*% cvfit$beta
      
      if (mcv == "Deviance"){
        
        test.tj <- sort(test.t[test.d==1])
        for (j in 1:length(test.tj)){
          cv.devnull <- cv.devnull + log(sum(test.t >= test.tj[j]))
          if (sum(test.t >= test.tj[j]) > 1){
            cvdev[,cv] <- cvdev[,cv] + linprod[test.t == test.tj[j],] - 
              log(colSums(exp(linprod[test.t >= test.tj[j],])) + 1e-8)
          }else if (sum(test.t >= test.tj[j]) == 1){
            cvdev[,cv] <- cvdev[,cv] + linprod[test.t == test.tj[j],] - 
              log(exp(linprod[test.t >= test.tj[j],]) + 1e-8) 
          }
          #- log(accu(link(widx))+1e-8)
        }
        
        cvdevnull[cv] <- 2*cv.devnull
        cvdev[,cv] <- -2*cvdev[,cv]
        
      }else if (mcv == "Cindex"){
        
        for (kk in 1:length(lidx)){
          cvdev[kk,cv] <- 1-concordance.index(x=linprod[,kk],surv.time=test.t,surv.event=test.d)$c.index
        }
        
      }
      
      
     
      # cvmse[,cv] <- apply(cbind(1,test.x) %*% rbind(t(cvfit$beta0),cvfit$beta),2,function(x) sum((test.y - exp(x)/(1+exp(x)))^2)/length(test.y))
      
    }
    
    mean.cvdev <- rowMeans(cvdev)
    se.cvdev <- apply(cvdev,1,sd)
    
    idx.min <- which.min(mean.cvdev)
    se.min <- se.cvdev[idx.min]
    idx.1se <- suppressWarnings(max(which(mean.cvdev > mean.cvdev[idx.min] + se.min & 1:length(lidx) < idx.min)))
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
      
      devplot <- data.frame(loglambda=log(ret$lambda),
                            dev=ret$cvdev.mean,
                            se=ret$cvdev.se,
                            devaddse=ret$cvdev.mean+ret$cvdev.se,
                            devminse=ret$cvdev.mean-ret$cvdev.se)
      
      pdev <- ggplot(devplot, aes(x=loglambda, y=dev)) +
        geom_errorbar(aes(ymin=devminse,ymax=devaddse),color="grey")+
        geom_point(color="red")+
        theme_bw() +
        xlab("log(lambda)") +
        ylab("Partial Likelihood Deviance")+
        geom_vline(xintercept=log(ret$lambda[ret$best.idx$idx.min]),linetype="dashed",color="darkgrey")+
        geom_vline(xintercept=log(ret$lambda[ret$best.idx$idx.1se]),linetype="dotted",color="darkgrey")+
        annotate("text",x=lambda_count$loglambda,y=max(devplot$devaddse)+0.05,label=as.character(lambda_count$count))+
        ggtitle("Cross-validated deviance versus log(lambda)")
      
      ret$pcoef <- pcoef
      ret$pdev <- pdev
      
    }
    
    if (step2){ # need to develop a equivalent lasso procedure for this. Stepwise selection is too slow for a big number of selected variables.
      
      if (length(which(ret$best.beta$min.mse!=0)) <= 10 & length(which(ret$best.beta$min.mse!=0)) > 0){
        idxs <- combn(which(ret$best.beta$min.mse!=0),2)
        x.select.min <- matrix(NA,nrow=n,ncol=ncol(idxs))
        for (k in 1:ncol(idxs)){
          x.select.min[,k] <- x[,idxs[1,k]] - x[,idxs[2,k]]
        }
        df_step2 <- data.frame(d=d,t=t,x=x.select.min)
        step2fit <- step(coxph(Surv(t,d)~.,data=df_step2),trace=0)
        vars <- as.numeric(sapply(names(step2fit$coefficients),function(x) strsplit(x,split = "[.]")[[1]][2]))
        
        if (ncol(idxs) == 1 & length(vars) == 2){
          vars = 1
        }
        
        selected <- idxs[,vars]
        for (k1 in 1:nrow(selected)){
          for (k2 in 1:ncol(selected)){
            selected[k1,k2] <- colnames(x)[as.numeric(selected[k1,k2])]
          }
        }
        ret$step2.feature.min = selected
        ret$step2fit.min <- step2fit
      }
      
      if (length(which(ret$best.beta$add.1se!=0)) <= 10 & length(which(ret$best.beta$add.1se!=0)) > 0){
        idxs <- combn(which(ret$best.beta$add.1se!=0),2)
        x.select.min <- matrix(NA,nrow=n,ncol=ncol(idxs))
        for (k in 1:ncol(idxs)){
          x.select.min[,k] <- x[,idxs[1,k]] - x[,idxs[2,k]]
        }
        df_step2 <- data.frame(d=d,t=t,x=x.select.min)
        step2fit <- step(coxph(Surv(t,d)~.,data=df_step2),trace=0)
        vars <- as.numeric(sapply(names(step2fit$coefficients),function(x) strsplit(x,split = "[.]")[[1]][2]))
        
        if (ncol(idxs) == 1 & length(vars) == 2){
          vars = 1
        }
        
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
    
    ret <- list(beta=fullfit$beta[,lidx],
                lambda=fullfit$lambda[lidx],
                loss=fullfit$loss[lidx],
                mse=fullfit$mse[lidx]
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