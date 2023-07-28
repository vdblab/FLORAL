LogRatioFGLasso <- function(x,
                            y,
                            id,
                            weight,
                            ncov,
                            length.lambda=100,
                            lambda.min.ratio=NULL,
                            wcov,
                            a=1,
                            mu=1,
                            ncv=5,
                            foldid=NULL,
                            step2=FALSE,
                            progress=TRUE,
                            plot=TRUE,
                            mcv="Deviance",
                            ncore=1){
  
  ptm <- proc.time()
  
  t0 <- y[,1]
  t1 <- y[,2]
  d <- y[,3]
  n <- nrow(y)
  p <- ncol(x)
  tt <- y[,1:2]
  
  tj <- sort(t1[d==1])
  if (length(unique(t1)) < length(t1)){
    t1[t1 %in% tj] <- t1[t1 %in% tj] + runif(length(t1[t1 %in% tj]),min=0,max=1e-4) # break ties
  }  
  tj <- sort(t1[d==1])
  
  devnull <- log(sum(weight*(t1 >= tj[1] & t0 < tj[1])))
  denomj <- rep(NA,length(tj))
  denomj[1] <- sum(weight*(t1 >= tj[1] & t0 < tj[1]))
  for (j in 2:length(tj)){
    denomj[j] <- sum(weight*(t1 >= tj[j] & t0 < tj[j]))
    devnull <- devnull + log(denomj[j])
  }
  devnull <- 2*devnull
  
  expect <- apply(tt,1,function(x) sum(1/denomj[which(tj <= x[2] & tj > x[1])]))
  sfun = d - expect
  # hfun <- expect - sapply(t,function(x) sum(1/denomj[which(tj <= x)]^2)) + 1e-8 # hessian
  # z <- sfun/hfun
  
  if (a > 0){
    lambda0 <- max(abs(t(sfun) %*% x))/(a*n)
  }else if (a == 0){
    lambda0 <- max(abs(t(sfun) %*% x))/(1e-3*n)
  }
  
  adjust = FALSE
  if (ncov > 0) adjust = TRUE
  
  if (is.null(lambda.min.ratio)) lambda.min.ratio = ifelse(n < p, 1e-02, 1e-02)
  lambda <- 10^(seq(log10(lambda0),log10(lambda0*lambda.min.ratio),length.out=length.lambda))
  
  if (progress) cat("Algorithm running for full dataset: \n")
  
  fullfit <- fg_enet_al(x,t0,t1,d,tj,weight,length.lambda,mu,100,lambda,wcov,a,adjust,ncov,devnull,progress)
  lidx <- which(fullfit$loglik != 0 | !is.nan(fullfit$loglik))
  
  dev <- -2*fullfit$loglik[lidx]
  
  if (!is.null(colnames(x))){
    rownames(fullfit$beta) = colnames(x)
  }else{
    colnames(x) = 1:ncol(x)
    rownames(fullfit$beta) = 1:ncol(x)
  }
  
  if (!is.null(ncv)){
    
    cvdev <- matrix(0,nrow=length(lidx),ncol=ncv)
    cvdevnull <- rep(0,ncol=ncv)
    
    if (is.null(foldid)){
      labels <- coxsplitss(as.matrix(y),id,ncv)
    }else{
      labels <- fgfoldid(id,foldid)
    }
    
    # labels <- caret::createFolds(factor(d),k=ncv)
    
    if (ncore == 1){
      
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
        train.t0 <- t0[labels!=cv]
        train.t1 <- t1[labels!=cv]
        train.tt <- y[labels!=cv,1:2]
        train.w <- weight[labels!=cv]
        
        test.x <- x[labels==cv,]
        test.d <- d[labels==cv]
        test.t0 <- t0[labels==cv]
        test.t1 <- t1[labels==cv]
        test.tt <- y[labels==cv,1:2]
        test.w <- weight[labels==cv]
        
        cv.devnull <- 0
        train.tj <- sort(train.t1[train.d==1])
        for (j in 1:length(train.tj)){
          cv.devnull <- cv.devnull + log(sum(train.t1 >= train.tj[1] & train.t0 < train.tj[1]))
        }
        cv.devnull <- 2*cv.devnull
        
        cvfit <- fg_enet_al(train.x,train.t0,train.t1,train.d,train.tj,train.w,length(lidx),mu,100,lambda[lidx],wcov,a,adjust,ncov,cv.devnull,progress)
        
        cv.devnull <- 0
        loglik <- rep(0,length(lidx))
        linprod <- test.w * (test.x %*% cvfit$beta)
        
        if (mcv == "Deviance"){
          
          test.tj <- sort(test.t1[test.d==1])
          for (j in 1:length(test.tj)){
            cv.devnull <- cv.devnull + log(sum(test.t1 >= test.tj[j] & test.t0 < test.tj[j]))
            if (sum(test.t1 >= test.tj[j] & test.t0 < test.tj[j]) > 1){
              cvdev[,cv] <- cvdev[,cv] + linprod[test.t1 == test.tj[j],] - 
                log(colSums(exp(linprod[test.t1 >= test.tj[j] & test.t0 < test.tj[j],])) + 1e-8)
            }else if (sum(test.t1 >= test.tj[j] & test.t0 < test.tj[j]) == 1){
              cvdev[,cv] <- cvdev[,cv] + linprod[test.t1 == test.tj[j],] - 
                log(exp(linprod[test.t1 >= test.tj[j] & test.t0 < test.tj[j],]) + 1e-8) 
            }
            #- log(accu(link(widx))+1e-8)
          }
          
          cvdevnull[cv] <- 2*cv.devnull
          cvdev[,cv] <- -2*cvdev[,cv]
          
        }
        
        # cvmse[,cv] <- apply(cbind(1,test.x) %*% rbind(t(cvfit$beta0),cvfit$beta),2,function(x) sum((test.y - exp(x)/(1+exp(x)))^2)/length(test.y))
        
      }
      
    }else if(ncore > 1){
      
      if (progress) warning(paste0("Using ", ncore ," core for cross-validation computation."))
      
      cl <- makeCluster(ncore)
      registerDoParallel(cl)
      
      cvdev <- foreach(cv=1:ncv,.combine=cbind) %dopar% {
        
        train.x <- x[labels!=cv,]
        train.d <- d[labels!=cv]
        train.t0 <- t0[labels!=cv]
        train.t1 <- t1[labels!=cv]
        train.tt <- y[labels!=cv,1:2]
        train.w <- weight[labels!=cv]
        
        test.x <- x[labels==cv,]
        test.d <- d[labels==cv]
        test.t0 <- t0[labels==cv]
        test.t1 <- t1[labels==cv]
        test.tt <- y[labels==cv,1:2]
        test.w <- weight[labels==cv]
        
        cv.devnull <- 0
        train.tj <- sort(train.t1[train.d==1])
        for (j in 1:length(train.tj)){
          cv.devnull <- cv.devnull + log(sum(train.t1 >= train.tj[1] & train.t0 < train.tj[1]))
        }
        cv.devnull <- 2*cv.devnull
        
        cvfit <- fg_enet_al(train.x,train.t0,train.t1,train.d,train.tj,train.w,length(lidx),mu,100,lambda[lidx],wcov,a,adjust,ncov,cv.devnull,FALSE)
        
        cv.dev <- rep(0,length(lidx))
        linprod <- test.w * (test.x %*% cvfit$beta)
        
        test.tj <- sort(test.t1[test.d==1])
        for (j in 1:length(test.tj)){
          if (sum(test.t1 >= test.tj[j] & test.t0 < test.tj[j]) > 1){
            cv.dev <- cv.dev + linprod[test.t1 == test.tj[j],] - 
              log(colSums(exp(linprod[test.t1 >= test.tj[j] & test.t0 < test.tj[j],])) + 1e-8)
          }else if (sum(test.t1 >= test.tj[j] & test.t0 < test.tj[j]) == 1){
            cv.dev <- cv.dev + linprod[test.t1 == test.tj[j],] - 
              log(exp(linprod[test.t1 >= test.tj[j] & test.t0 < test.tj[j],]) + 1e-8) 
          }
        }
        
        cv.dev <- -2*cv.dev
        cv.dev
        
      }
      
      stopCluster(cl)
      
    }
    
    mean.cvdev <- rowMeans(cvdev)
    se.cvdev <- apply(cvdev,1,function(x) sd(x)/sqrt(ncv))
    
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
                a=a,
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
      
      devplot <- data.frame(loglambda=log(ret$lambda),
                            dev=ret$cvdev.mean,
                            se=ret$cvdev.se,
                            devaddse=ret$cvdev.mean+ret$cvdev.se,
                            devminse=ret$cvdev.mean-ret$cvdev.se)
      
      pdev <- ggplot(devplot, aes(x=.data$loglambda, y=.data$dev)) +
        geom_errorbar(aes(ymin=.data$devminse,ymax=.data$devaddse),color="grey")+
        geom_point(color="red")+
        theme_bw() +
        xlab(expression(paste("log(",lambda,")"))) +
        ylab("Partial Likelihood Deviance")+
        geom_vline(xintercept=log(ret$lambda[ret$best.idx$idx.min]),linetype="dashed",color="darkgrey")+
        geom_vline(xintercept=log(ret$lambda[ret$best.idx$idx.1se]),linetype="dotted",color="darkgrey")+
        annotate("text",x=lambda_count$loglambda,y=max(devplot$devaddse)+0.05,label=as.character(lambda_count$count))+
        ggtitle(expression(paste("Cross-validated MSE versus log(",lambda,")")))
      
      ret$pcoef <- pcoef
      ret$pdev <- pdev
      
    }
    
    if (step2){ # need to develop a equivalent lasso procedure for this. Stepwise selection is too slow for a big number of selected variables.
      
      if (!adjust){
        
        if(progress) cat("Step 2 filtering started.\n")
        
        if (length(which(ret$best.beta$min.mse!=0)) > 1){
          
          idxs <- combn(which(ret$best.beta$min.mse!=0),2)
          
          x.select.min <- matrix(NA,nrow=n,ncol=ncol(idxs))
          for (k in 1:ncol(idxs)){
            x.select.min[,k] <- x[,idxs[1,k]] - x[,idxs[2,k]]
          }
          
          if (ncol(x.select.min) > 1){
            stepglmnet <- suppressWarnings(cv.glmnet(x=x.select.min,y=Surv(t0,t1,d),weights=weight,type.measure = "deviance",family="cox"))
            x.select.min <- x.select.min[,which(stepglmnet$glmnet.fit$beta[,stepglmnet$index[1]]!=0)]
            idxs <- idxs[,which(stepglmnet$glmnet.fit$beta[,stepglmnet$index[1]]!=0)]
          }else{
            idxs <- as.vector(idxs)
          }
          
          df_step2 <- data.frame(t0=t0,t1=t1,d=d,x=x.select.min)
          step2fit <- suppressWarnings(step(coxph(Surv(t0,t1,d)~.,weights=weight,data=df_step2),trace=0))
          vars <- as.numeric(sapply(names(step2fit$coefficients),function(x) strsplit(x,split = "[.]")[[1]][2]))
          
          if (is.null(ncol(idxs))){
            selected <- idxs
          }else{
            selected <- idxs[,vars]
          }
          # for (k1 in 1:nrow(selected)){
          #   for (k2 in 1:ncol(selected)){
          #     selected[k1,k2] <- colnames(x)[as.numeric(selected[k1,k2])]
          #   }
          # }
          
          ret$step2.feature.min = selected
          ret$step2fit.min <- step2fit
        }
        
        if (length(which(ret$best.beta$add.1se!=0)) > 1){
          
          idxs <- combn(which(ret$best.beta$add.1se!=0),2)
          
          x.select.min <- matrix(NA,nrow=n,ncol=ncol(idxs))
          for (k in 1:ncol(idxs)){
            x.select.min[,k] <- x[,idxs[1,k]] - x[,idxs[2,k]]
          }
          
          if (ncol(x.select.min) > 1){
            stepglmnet <- suppressWarnings(cv.glmnet(x=x.select.min,y=Surv(t0,t1,d),weights=weight,type.measure = "deviance",family="cox"))
            x.select.min <- x.select.min[,which(stepglmnet$glmnet.fit$beta[,stepglmnet$index[1]]!=0)]
            idxs <- idxs[,which(stepglmnet$glmnet.fit$beta[,stepglmnet$index[1]]!=0)]
          }else{
            idxs <- as.vector(idxs)
          }
          
          df_step2 <- data.frame(t0=t0,t1=t1,d=d,x=x.select.min)
          step2fit <- suppressWarnings(step(coxph(Surv(t0,t1,d)~.,weights=weight,data=df_step2),trace=0))
          vars <- as.numeric(sapply(names(step2fit$coefficients),function(x) strsplit(x,split = "[.]")[[1]][2]))
          
          if (is.null(ncol(idxs))){
            selected <- idxs
          }else{
            selected <- idxs[,vars]
          }
          # for (k1 in 1:nrow(selected)){
          #   for (k2 in 1:ncol(selected)){
          #     selected[k1,k2] <- colnames(x)[as.numeric(selected[k1,k2])]
          #   }
          # }
          
          ret$step2.feature.1se = selected
          ret$step2fit.1se <- step2fit
        }
        
      }else{
        
        if (length(which(ret$best.beta$min.mse!=0)) > 1){
          
          allidx <- which(ret$best.beta$min.mse!=0)
          
          covidx <- allidx[allidx <= ncov]
          taxidx <- allidx[allidx > ncov]
          
          idxs <- NULL
          x.select.min <- NULL
          
          if (length(taxidx) > 1){
            idxs <- combn(taxidx,2)
            for (k in 1:ncol(idxs)){
              x.select.min <- cbind(x.select.min, x[,idxs[1,k]] - x[,idxs[2,k]])
            }
            colnames(x.select.min) <- rep("",ncol(x.select.min))
          }
          
          if (length(covidx) > 0){
            x.select.min <- cbind(x.select.min, x[,covidx])
            if (!is.null(idxs)){
              colnames(x.select.min)[(ncol(idxs)+1):ncol(x.select.min)] = colnames(x)[covidx]
            }else{
              colnames(x.select.min) <- colnames(x)[covidx]
            }
          }
          
          
          # if(is.null(x.select.min)) break
          
          if (ncol(x.select.min) > 1){
            stepglmnet <- cv.glmnet(x=x.select.min,y=Surv(t0,t1,d),weights=weight,type.measure = "deviance",family="cox")
            
            if (length(taxidx) < 2){
              
              idxs <- NULL
              
              # covs <- colnames(x.select.min)[which(stepglmnet$glmnet.fit$beta[,stepglmnet$index[1]]!=0)]
              
            }else{
              
              if (length(covidx) == 0)  idxs <- idxs[,which(stepglmnet$glmnet.fit$beta[,stepglmnet$index[1]]!=0)]
              
              if (length(covidx) > 0){
                
                # covs <- colnames(x.select.min)[setdiff(which(stepglmnet$glmnet.fit$beta[,stepglmnet$index[1]]!=0),
                #                                        1:ncol(idxs))]
                
                idxs <- idxs[,setdiff(which(stepglmnet$glmnet.fit$beta[,stepglmnet$index[1]]!=0),
                                      (ncol(idxs)+1):(ncol(idxs)+length(covidx)))]
                
              }
              
            }
            
            x.select.min <- x.select.min[,which(stepglmnet$glmnet.fit$beta[,stepglmnet$index[1]]!=0)]
          }
          
          if (sum(colnames(x.select.min)=="") > 0)  colnames(x.select.min)[colnames(x.select.min)==""] <- paste0("x.",1:sum(colnames(x.select.min)==""))
          df_step2 <- data.frame(t0=t0,t1=t1,d=d,x.select.min)
          step2fit <- suppressWarnings(step(coxph(Surv(t0,t1,d)~.,weights=weight,data=df_step2),trace=0))
          vars <- suppressWarnings(as.numeric(sapply(names(step2fit$coefficients),function(x) strsplit(x,split = "[.]")[[1]][2])))
          vars <- vars[!is.na(vars)]
          
          selected <- NULL
          if (is.null(ncol(idxs)) & length(vars) > 0){
            selected <- idxs
          }else if (length(vars) > 0){
            selected <- idxs[,vars]
          }
          
          ret$step2.feature.min = selected
          ret$step2fit.min <- step2fit
        }
        
        if (length(which(ret$best.beta$add.1se!=0)) > 1){
          
          allidx <- which(ret$best.beta$add.1se!=0)
          
          covidx <- allidx[allidx <= ncov]
          taxidx <- allidx[allidx > ncov]
          
          idxs <- NULL
          x.select.min <- NULL
          
          if (length(taxidx) > 1){
            idxs <- combn(taxidx,2)
            for (k in 1:ncol(idxs)){
              x.select.min <- cbind(x.select.min, x[,idxs[1,k]] - x[,idxs[2,k]])
            }
            colnames(x.select.min) <- rep("",ncol(x.select.min))
          }
          
          if (length(covidx) > 0){
            x.select.min <- cbind(x.select.min, x[,covidx])
            if (!is.null(idxs)){
              colnames(x.select.min)[(ncol(idxs)+1):ncol(x.select.min)] = colnames(x)[covidx]
            }else{
              colnames(x.select.min) <- colnames(x)[covidx]
            }
          }
          
          # if(is.null(x.select.min)) break
          
          if (ncol(x.select.min) > 1){
            stepglmnet <- cv.glmnet(x=x.select.min,y=Surv(t0,t1,d),weights=weight,type.measure = "deviance",family="cox")
            
            if (length(taxidx) < 2){
              
              idxs <- NULL
              
              # covs <- colnames(x.select.min)[which(stepglmnet$glmnet.fit$beta[,stepglmnet$index[1]]!=0)]
              
            }else{
              
              if (length(covidx) == 0)  idxs <- idxs[,which(stepglmnet$glmnet.fit$beta[,stepglmnet$index[1]]!=0)]
              
              if (length(covidx) > 0){
                
                # covs <- colnames(x.select.min)[setdiff(which(stepglmnet$glmnet.fit$beta[,stepglmnet$index[1]]!=0),1:ncol(idxs))]
                
                idxs <- idxs[,setdiff(which(stepglmnet$glmnet.fit$beta[,stepglmnet$index[1]]!=0),
                                      (ncol(idxs)+1):(ncol(idxs)+length(covidx)))]
                
              }
              
            }
            
            x.select.min <- x.select.min[,which(stepglmnet$glmnet.fit$beta[,stepglmnet$index[1]]!=0)]
          }
          
          if (sum(colnames(x.select.min)=="") > 0)  colnames(x.select.min)[colnames(x.select.min)==""] <- paste0("x.",1:sum(colnames(x.select.min)==""))
          df_step2 <- data.frame(t0=t0,t1=t1,d=d,x.select.min)
          step2fit <- suppressWarnings(step(coxph(Surv(t0,t1,d)~.,weights=weight,data=df_step2),trace=0))
          vars <- suppressWarnings(as.numeric(sapply(names(step2fit$coefficients),function(x) strsplit(x,split = "[.]")[[1]][2])))
          vars <- vars[!is.na(vars)]
          
          selected <- NULL
          if (is.null(ncol(idxs)) & length(vars) > 0){
            selected <- idxs
          }else if (length(vars) > 0){
            selected <- idxs[,vars]
          }
          
          ret$step2.feature.1se = selected
          ret$step2fit.1se <- step2fit
          
        }
        
      }
      
    }
    
  }else{
    
    ret <- list(beta=fullfit$beta[,lidx],
                lambda=fullfit$lambda[lidx],
                a=a,
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