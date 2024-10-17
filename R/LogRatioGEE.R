LogRatioGEE <- function(x,
                        y,
                        id,
                        ncov, #not yet + intercept option (seems fine now)
                        intercept,
                        family,
                        corstr,
                        scalefix,
                        scalevalue,
                        length.lambda=100,
                        lambda.min.ratio=NULL,
                        wcov, # not yet (seems fine now)
                        a=1, # not yet finished development (seems fine now)
                        mu=1e6,
                        ncv=5,
                        foldid=NULL,
                        step2=FALSE, # not yet
                        progress=TRUE,
                        plot=TRUE,
                        ncore=1){
  
  ptm <- proc.time()
  
  n <- length(unique(id))
  p <- ncol(x)
  
  if (family == "gaussian"){
    
    if (a > 0){
      lambda0 <- max(abs(t(scale(y)) %*% x))/(min(a,1)*nrow(x))
    }else if (a == 0){
      lambda0 <- max(abs(t(scale(y)) %*% x))/(1e-3*nrow(x))
    }
    
  }else if (family == "binomial"){
    
    sfun = y-0.5
    
    if (a > 0){
      lambda0 <- max(abs(t(sfun) %*% x))/(min(a,1)*nrow(x))
    }else if (a == 0){
      lambda0 <- max(abs(t(sfun) %*% x))/(1e-3*nrow(x))
    }
    
  }
  
  family0 <- family
  if (is.character(family)) family <- get(family)
  if (is.function(family))  family <- family()
  
  if (is.null(lambda.min.ratio)) lambda.min.ratio = ifelse(n < p, 1e-02, 1e-02)
  lambda <- 10^(seq(log10(lambda0),log10(lambda0*lambda.min.ratio),length.out=length.lambda))
  
  nt <- as.integer(unlist(lapply(split(id, id), "length"))) # Number of obs per subject
  
  if (progress) cat("Algorithm running for full dataset: \n")
  
  fullfit <- gee_fit(y,
                     x,
                     nt, 
                     family$linkinv,
                     family$mu.eta,
                     family$variance,
                     corstr,
                     lambda,
                     a,
                     ncov,
                     wcov,
                     tol=1e-3,
                     eps=1e-6,
                     muu=mu,
                     maxiter1=100,
                     maxiter2=1,
                     scalefix=scalefix,
                     scalevalue=scalevalue,
                     display_progress=progress)
  
  if (!is.null(colnames(x))){
    rownames(fullfit$beta) = colnames(x)
  }else{
    colnames(x) = 1:ncol(x)
    rownames(fullfit$beta) = 1:ncol(x)
  }
  
  beta_filtered <- fullfit$beta
  beta_filtered[abs(beta_filtered) < 1e-3] = 0
  beta_filtered[apply(beta_filtered, 2,function(x) abs(x) < max(abs(x))*0.01)] <- 0
  
  if (!is.null(ncv)){
    
    cvmse <- matrix(NA,nrow=length.lambda,ncol=ncv)
    
    if (is.null(foldid)){
      
      ### Double check the way to split folds!
      
      labels <- caret::createFolds(unique(id),k=ncv,list=FALSE)
    }else{
      labels <- foldid
    }
    
    if (ncore == 1){
      
      for (cv in 1:ncv){
        
        if (progress) cat(paste0("Algorithm running for cv dataset ",cv," out of ",ncv,": \n"))
        
        train.x <- x[labels[id]!=cv,]
        train.y <- y[labels[id]!=cv]
        test.x <- x[labels[id]==cv,]
        test.y <- y[labels[id]==cv]
        train.nt <- nt[labels!=cv]
        test.nt <- nt[labels==cv]
        
        cvfit <- gee_fit(train.y,
                         train.x,
                         train.nt, 
                         family$linkinv,
                         family$mu.eta,
                         family$variance,
                         corstr,
                         lambda,
                         a,
                         ncov,
                         wcov,
                         tol=1e-3,
                         eps=1e-6,
                         muu=mu,
                         maxiter1=100,
                         maxiter2=1,
                         scalefix=scalefix,
                         scalevalue=scalevalue,
                         display_progress=progress)
        
        mufit=family$linkinv(test.x %*% cvfit$beta)
        
        #### !check
        
        cvmse[,cv] <- apply(mufit,2,function(xx) sum(family$dev.resids(test.y,xx,wt=1)))
        
      }
      
    }else if(ncore > 1){
      
      if (progress) cat(paste0("Using ", ncore ," core for cross-validation computation.\n"))
      
      Sys.sleep(1)
      
      cl <- makeCluster(ncore)
      registerDoParallel(cl)
      
      cvmse <- foreach(cv=1:ncv,.combine=cbind) %dopar% {
        
        train.x <- x[labels[id]!=cv,]
        train.y <- y[labels[id]!=cv]
        test.x <- x[labels[id]==cv,]
        test.y <- y[labels[id]==cv]
        train.nt <- nt[labels!=cv]
        test.nt <- nt[labels==cv]
        
        cvfit <- gee_fit(train.y,
                         train.x,
                         train.nt, 
                         family$linkinv,
                         family$mu.eta,
                         family$variance,
                         corstr,
                         lambda,
                         a,
                         ncov,
                         wcov,
                         tol=1e-3,
                         eps=1e-6,
                         muu=mu,
                         maxiter1=100,
                         maxiter2=1,
                         scalefix=scalefix,
                         scalevalue=scalevalue,
                         display_progress=progress)
        
        mufit=family$linkinv(test.x %*% cvfit$beta)
        apply(mufit,2,function(x) sum(family$dev.resids(test.y,x,wt=1)))
        
      }
      
      stopCluster(cl)
      
    }
    
    mean.cvmse <- rowMeans(cvmse)
    se.cvmse <- apply(cvmse,1,function(x) sd(x)/sqrt(ncv))
    median.cvmse <- apply(cvmse,1,median)
    
    idx.min <- which.min(mean.cvmse)
    se.min <- se.cvmse[idx.min]
    idx.1se <- suppressWarnings(min(which(mean.cvmse < mean.cvmse[idx.min] + se.min & 1:length.lambda < idx.min)))
    if (idx.1se == Inf) idx.1se = idx.min
    
    best.beta <- list(min.mse = beta_filtered[,idx.min],
                      add.1se = beta_filtered[,idx.1se])
    
    best.idx <- list(idx.min = idx.min,
                     idx.1se = idx.1se)
    
    ret <- list(beta=beta_filtered,
                beta_unfiltered = fullfit$beta,
                lambda=lambda,
                a=a,
                cvmse.mean=mean.cvmse,
                cvmse.se=se.cvmse,
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
      
      pcoef <- ggplot(beta_nzero, aes(x=.data$loglambda,y=.data$value,group=.data$X1,color=as.factor(.data$X1))) + 
        geom_line() + 
        scale_color_manual(values=rainbow(sum(rowSums(ret$beta != 0) > 0))) + 
        theme_bw() + 
        theme(legend.position = "none") +
        xlab(expression(paste("log(",lambda,")"))) +
        ylab("Coefficient") +
        geom_vline(xintercept=log(ret$lambda[ret$best.idx$idx.min]),linetype="dashed",color="darkgrey")+
        geom_vline(xintercept=log(ret$lambda[ret$best.idx$idx.1se]),linetype="dotted",color="darkgrey")+
        # annotate("text",x=min(beta_nzero$loglambda)-2,y=top10feat,label=top10name,hjust=0)+
        annotate("text",x=lambda_count$loglambda,y=max(beta_nzero$value)+0.2,label=as.character(lambda_count$count))+
        ggtitle(expression(paste("Coefficients versus log(",lambda,")")))
      
      mseplot <- data.frame(loglambda=log(ret$lambda),
                            mse=ret$cvmse.mean,
                            se=ret$cvmse.se,
                            mseaddse=ret$cvmse.mean+ret$cvmse.se,
                            mseminse=ret$cvmse.mean-ret$cvmse.se)
      
      pmse <- ggplot(mseplot, aes(x=.data$loglambda, y=.data$mse)) +
        geom_errorbar(aes(ymin=.data$mseminse,ymax=.data$mseaddse),color="grey")+
        geom_point(color="red")+
        theme_bw() +
        xlab(expression(paste("log(",lambda,")"))) +
        ylab("Mean-Squared Error")+
        geom_vline(xintercept=log(ret$lambda[ret$best.idx$idx.min]),linetype="dashed",color="darkgrey")+
        geom_vline(xintercept=log(ret$lambda[ret$best.idx$idx.1se]),linetype="dotted",color="darkgrey")+
        annotate("text",x=lambda_count$loglambda,y=max(mseplot$mseaddse)+0.05,label=as.character(lambda_count$count))+
        ggtitle(expression(paste("Cross-validated deviance residual versus log(",lambda,")")))
      
      ret$pcoef <- pcoef
      ret$pmse <- pmse
      
    }
    
    if (step2){
      
      if (progress) cat(paste0("Step2 started."))
      
      if (ncov > 0){
        idxfeat <- setdiff(1:length(ret$best.beta$min.mse),1:ncov)
      }else{
        idxfeat <- 1:length(ret$best.beta$min.mse)
      }
      
      if (length(which(ret$best.beta$min.mse[idxfeat]!=0)) > 1){
        # idxs <- combn(which(ret$best.beta$min.mse[idxfeat]!=0),2)
        idxs <-combn(names(which(ret$best.beta$min.mse[idxfeat]!=0)),2)
        
        x.select.min <- matrix(NA,nrow=nrow(x),ncol=ncol(idxs))
        for (k in 1:ncol(idxs)){
          x.select.min[,k] <- x[,idxs[1,k]] - x[,idxs[2,k]]
        }
        
        if (ncov > 0){
          x.select.min <- cbind(x[,1:ncov],x.select.min)
        }
        
        if (ncol(x.select.min) > ncov+1){
          
          if (family0 == "gaussian"){
            
            if (a > 0){
              lambda0 <- max(abs(t(scale(y)) %*% x.select.min))/(min(a,1)*nrow(x.select.min))
            }else if (a == 0){
              lambda0 <- max(abs(t(scale(y)) %*% x.select.min))/(1e-3*nrow(x.select.min))
            }
            
          }else if (family0 == "binomial"){
            
            sfun = y-0.5
            
            if (a > 0){
              lambda0 <- max(abs(t(sfun) %*% x.select.min))/(min(a,1)*nrow(x.select.min))
            }else if (a == 0){
              lambda0 <- max(abs(t(sfun) %*% x.select.min))/(1e-3*nrow(x.select.min))
            }
            
          }
          
          # if (is.null(lambda.min.ratio)) 
          lambda.min.ratio = ifelse(n < p, 1e-02, 1e-02)
          lambda <- 10^(seq(log10(lambda0),log10(lambda0*lambda.min.ratio),length.out=length.lambda))
          
          fullfit <- gee_fit(y,
                             x.select.min,
                             nt, 
                             family$linkinv,
                             family$mu.eta,
                             family$variance,
                             corstr,
                             lambda,
                             a,
                             ncov,
                             wcov,
                             tol=1e-3,
                             eps=1e-6,
                             muu=0,
                             maxiter1=100,
                             maxiter2=1,
                             scalefix=scalefix,
                             scalevalue=scalevalue,
                             display_progress=FALSE)
          
          cvmse <- matrix(NA,nrow=length.lambda,ncol=ncv)
          
          if (ncore == 1){
            
            for (cv in 1:ncv){
              
              train.x <- x.select.min[labels[id]!=cv,]
              train.y <- y[labels[id]!=cv]
              test.x <- x.select.min[labels[id]==cv,]
              test.y <- y[labels[id]==cv]
              train.nt <- nt[labels!=cv]
              test.nt <- nt[labels==cv]
              
              cvfit <- gee_fit(train.y,
                               train.x,
                               train.nt, 
                               family$linkinv,
                               family$mu.eta,
                               family$variance,
                               corstr,
                               lambda,
                               a,
                               ncov,
                               wcov,
                               tol=1e-3,
                               eps=1e-6,
                               muu=0,
                               maxiter1=100,
                               maxiter2=1,
                               scalefix=scalefix,
                               scalevalue=scalevalue,
                               display_progress=FALSE)
              
              mufit=family$linkinv(test.x %*% cvfit$beta)
              cvmse[,cv] <- apply(mufit,2,function(x) sum(family$dev.resids(test.y,x,wt=1)))
              
            }
            
          }else if(ncore > 1){
            
            if (progress) cat(paste0("Step2: Using ", ncore ," core for cross-validation computation."))
            
            Sys.sleep(1)
            
            cl <- makeCluster(ncore)
            registerDoParallel(cl)
            
            cvmse <- foreach(cv=1:ncv,.combine=cbind) %dopar% {
              
              train.x <- x.select.min[labels[id]!=cv,]
              train.y <- y[labels[id]!=cv]
              test.x <- x.select.min[labels[id]==cv,]
              test.y <- y[labels[id]==cv]
              train.nt <- nt[labels!=cv]
              test.nt <- nt[labels==cv]
              
              cvfit <- gee_fit(train.y,
                               train.x,
                               train.nt, 
                               family$linkinv,
                               family$mu.eta,
                               family$variance,
                               corstr,
                               lambda,
                               a,
                               ncov,
                               wcov,
                               tol=1e-3,
                               eps=1e-6,
                               muu=0,
                               maxiter1=100,                         
                               maxiter2=1,
                               scalefix=scalefix,
                               scalevalue=scalevalue,
                               display_progress=FALSE)
              
              mufit=family$linkinv(test.x %*% cvfit$beta)
              apply(mufit,2,function(x) sum(family$dev.resids(test.y,x,wt=1)))
              
            }
            
            stopCluster(cl)
            
          }
          
          mean.cvmse <- rowMeans(cvmse)
          se.cvmse <- apply(cvmse,1,function(x) sd(x)/sqrt(ncv))
          
          idx.min <- which.min(mean.cvmse)
          se.min <- se.cvmse[idx.min]
          idx.1se <- suppressWarnings(min(which(mean.cvmse < mean.cvmse[idx.min] + se.min & 1:length.lambda < idx.min)))
          
          betafilt <- fullfit$beta
          betafilt[abs(betafilt) < 1e-3] = 0
          betafilt[apply(betafilt, 2,function(x) abs(x) < max(abs(x))*0.01)] <- 0
          
          # x.select.min <- x.select.min[,which(betafilt[,idx.1se]!=0)]
          
          taxanames <- as.vector(na.omit(apply(idxs,2,function(x) ifelse(sum(is.na(x))==0,paste(x,collapse ="/"),NA))))
          
          if (idx.1se == Inf){
            if (ncov > 0){
              idxs <- idxs[,which(betafilt[setdiff(1:length(betafilt[,idx.min]),1:ncov),idx.min]!=0)]
            }else{
              idxs <- idxs[,which(betafilt[,idx.min]!=0)]
            }
            
            coefs <- betafilt[,idx.min]
            if (ncov > 0){
              names(coefs)[1:ncov] <- colnames(x)[1:ncov]
            }
            names(coefs)[(ncov+1):(ncov+length(taxanames))] <- taxanames
            coefs <- coefs[coefs!=0]
            
          }else{
            if (ncov > 0){
              idxs <- idxs[,which(betafilt[setdiff(1:length(betafilt[,idx.1se]),1:ncov),idx.1se]!=0)]
            }else{
              idxs <- idxs[,which(betafilt[,idx.1se]!=0)]
            }
            
            coefs <- betafilt[,idx.1se]
            if (ncov > 0){
              names(coefs)[1:ncov] <- colnames(x)[1:ncov]
            }
            names(coefs)[(ncov+1):(ncov+length(taxanames))] <- taxanames
            coefs <- coefs[coefs!=0]
          }
          
        }else{
          
          idxs <- as.vector(idxs)
          
          coefs <- NA
          
        }
        
        ret$step2.feature.min = idxs
        ret$step2fit.min <- coefs
        
      }
      
      if (ncov > 0){
        idxfeat <- setdiff(1:length(ret$best.beta$add.1se),1:ncov)
      }else{
        idxfeat <- 1:length(ret$best.beta$add.1se)
      }
      
      if (length(which(ret$best.beta$add.1se[idxfeat]!=0)) > 1){
        
        idxs <-combn(names(which(ret$best.beta$add.1se[idxfeat]!=0)),2)
        
        x.select.min <- matrix(NA,nrow=nrow(x),ncol=ncol(idxs))
        for (k in 1:ncol(idxs)){
          x.select.min[,k] <- x[,idxs[1,k]] - x[,idxs[2,k]]
        }
        
        if (ncov > 0){
          x.select.min <- cbind(x[,1:ncov],x.select.min)
        }
        
        if (ncol(x.select.min) > ncov+1){
          
          if (family0 == "gaussian"){
            
            if (a > 0){
              lambda0 <- max(abs(t(scale(y)) %*% x.select.min))/(min(a,1)*nrow(x.select.min))
            }else if (a == 0){
              lambda0 <- max(abs(t(scale(y)) %*% x.select.min))/(1e-3*nrow(x.select.min))
            }
            
          }else if (family0 == "binomial"){
            
            sfun = y-0.5
            
            if (a > 0){
              lambda0 <- max(abs(t(sfun) %*% x.select.min))/(min(a,1)*nrow(x.select.min))
            }else if (a == 0){
              lambda0 <- max(abs(t(sfun) %*% x.select.min))/(1e-3*nrow(x.select.min))
            }
            
          }
          
          # if (is.null(lambda.min.ratio)) 
          lambda.min.ratio = ifelse(n < p, 1e-02, 1e-02)
          lambda <- 10^(seq(log10(lambda0),log10(lambda0*lambda.min.ratio),length.out=length.lambda))
          
          fullfit <- gee_fit(y,
                             x.select.min,
                             nt, 
                             family$linkinv,
                             family$mu.eta,
                             family$variance,
                             corstr,
                             lambda,
                             a,
                             ncov,
                             wcov,
                             tol=1e-3,
                             eps=1e-6,
                             muu=0,
                             maxiter1=100,                          
                             maxiter2=1,
                             scalefix=scalefix,
                             scalevalue=scalevalue,
                             display_progress=FALSE)
          
          cvmse <- matrix(NA,nrow=length.lambda,ncol=ncv)
          
          if (ncore == 1){
            
            for (cv in 1:ncv){
              
              # if (progress) cat(paste0("Step2: Algorithm running for cv dataset ",cv," out of ",ncv,": \n"))
              
              train.x <- x.select.min[labels[id]!=cv,]
              train.y <- y[labels[id]!=cv]
              test.x <- x.select.min[labels[id]==cv,]
              test.y <- y[labels[id]==cv]
              train.nt <- nt[labels!=cv]
              test.nt <- nt[labels==cv]
              
              cvfit <- gee_fit(train.y,
                               train.x,
                               train.nt, 
                               family$linkinv,
                               family$mu.eta,
                               family$variance,
                               corstr,
                               lambda,
                               a,
                               ncov,
                               wcov,
                               tol=1e-3,
                               eps=1e-6,
                               muu=0,
                               maxiter1=100,                          
                               maxiter2=1,
                               scalefix=scalefix,
                               scalevalue=scalevalue,
                               display_progress=FALSE)
              
              mufit=family$linkinv(test.x %*% cvfit$beta)
              cvmse[,cv] <- apply(mufit,2,function(x) sum(family$dev.resids(test.y,x,wt=1)))
              
            }
            
          }else if(ncore > 1){
            
            if (progress) cat(paste0("Step2: Using ", ncore ," core for cross-validation computation."))
            
            Sys.sleep(1)
            
            cl <- makeCluster(ncore)
            registerDoParallel(cl)
            
            cvmse <- foreach(cv=1:ncv,.combine=cbind) %dopar% {
              
              train.x <- x.select.min[labels[id]!=cv,]
              train.y <- y[labels[id]!=cv]
              test.x <- x.select.min[labels[id]==cv,]
              test.y <- y[labels[id]==cv]
              train.nt <- nt[labels!=cv]
              test.nt <- nt[labels==cv]
              
              cvfit <- gee_fit(train.y,
                               train.x,
                               train.nt, 
                               family$linkinv,
                               family$mu.eta,
                               family$variance,
                               corstr,
                               lambda,
                               a,
                               ncov,
                               wcov,
                               tol=1e-3,
                               eps=1e-6,
                               muu=0,
                               maxiter1=100,                         
                               maxiter2=1,
                               scalefix=scalefix,
                               scalevalue=scalevalue,
                               display_progress=FALSE)
              
              mufit=family$linkinv(test.x %*% cvfit$beta)
              apply(mufit,2,function(x) sum(family$dev.resids(test.y,x,wt=1)))
              
            }
            
            stopCluster(cl)
            
          }
          
          mean.cvmse <- rowMeans(cvmse)
          se.cvmse <- apply(cvmse,1,function(x) sd(x)/sqrt(ncv))
          
          idx.min <- which.min(mean.cvmse)
          se.min <- se.cvmse[idx.min]
          idx.1se <- suppressWarnings(min(which(mean.cvmse < mean.cvmse[idx.min] + se.min & 1:length.lambda < idx.min)))

          betafilt <- fullfit$beta
          # betafilt[abs(betafilt) < 5e-3] = 0
          beta_filtered[abs(beta_filtered) < 1e-3] = 0
          beta_filtered[apply(beta_filtered, 2,function(x) abs(x) < max(abs(x))*0.01)] <- 0
          
          # x.select.min <- x.select.min[,which(betafilt[,idx.1se]!=0)]
          # idxs <- idxs[,which(betafilt[,idx.1se]!=0)]
          
          taxanames <- as.vector(na.omit(apply(idxs,2,function(x) ifelse(sum(is.na(x))==0,paste(x,collapse ="/"),NA))))
          
          if (idx.1se == Inf){
            if (ncov > 0){
              idxs <- idxs[,which(betafilt[setdiff(1:length(betafilt[,idx.min]),1:ncov),idx.min]!=0)]
            }else{
              idxs <- idxs[,which(betafilt[,idx.min]!=0)]
            }
            
            coefs <- betafilt[,idx.min]
            if (ncov > 0){
              names(coefs)[1:ncov] <- colnames(x)[1:ncov]
            }
            names(coefs)[(ncov+1):(ncov+length(taxanames))] <- taxanames
            coefs <- coefs[coefs!=0]
            
          }else{
            if (ncov > 0){
              idxs <- idxs[,which(betafilt[setdiff(1:length(betafilt[,idx.1se]),1:ncov),idx.1se]!=0)]
            }else{
              idxs <- idxs[,which(betafilt[,idx.1se]!=0)]
            }
            
            coefs <- betafilt[,idx.1se]
            if (ncov > 0){
              names(coefs)[1:ncov] <- colnames(x)[1:ncov]
            }
            names(coefs)[(ncov+1):(ncov+length(taxanames))] <- taxanames
            coefs <- coefs[coefs!=0]
          }
          
        }else{
          idxs <- as.vector(idxs)
          
          coefs <- NA
        }
        
        ret$step2.feature.1se = idxs
        ret$step2fit.1se <- coefs
        
      }
      
      if (progress) cat("Step2 is completed.")
    }
  }else{
    
    ret <- list(beta=beta_filtered,
                beta_unfiltered = fullfit$beta,
                lambda=lambda,
                a=a
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
      
      pcoef <- ggplot(beta_nzero, aes(x=.data$loglambda,y=.data$value,group=.data$X1,color=as.factor(.data$X1))) + 
        geom_line() + 
        scale_color_manual(values=rainbow(sum(rowSums(ret$beta != 0) > 0))) + 
        theme_bw() + 
        theme(legend.position = "none") +
        xlab(expression(paste("log(",lambda,")"))) +
        ylab("Coefficient") +
        # annotate("text",x=min(beta_nzero$loglambda)-2,y=top10feat,label=top10name,hjust=0)+
        annotate("text",x=lambda_count$loglambda,y=max(beta_nzero$value)+0.2,label=as.character(lambda_count$count))+
        ggtitle(expression(paste("Coefficients versus log(",lambda,")")))
      
      ret$pcoef <- pcoef
      
    }
    
  }
  
  ret$runtime <- proc.time() - ptm
  
  return(ret)
  
}