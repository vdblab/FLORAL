LogRatioLogisticLassoKN <- function(x,
                                    y,
                                    ncov,
                                    length.lambda=100,
                                    lambda.min.ratio=NULL,
                                    wcov,
                                    a=1,
                                    mu=1,
                                    fdr=0.05,
                                    kn_method,
                                    offset=0,
                                    maxiter=100,
                                    intercept=FALSE,
                                    progress=TRUE,
                                    plot=TRUE
                                    ){
  
  ptm <- proc.time()
  
  n <- length(y)
  sfun = y-0.5
  
  adjust = FALSE
  if (ncov > 0) adjust = TRUE
  
  if (adjust){
    
    x0 <- x[,(ncov+1):ncol(x)]
    
    if (kn_method$model == "2order"){
      xk <- knockoff::create.second_order(X=x0)
    }else if (kn_method$model == "VAE"){
      xk <- VAE_func(x0)$knockoff_x
    }else if (kn_method$model == "DK"){
      xk <- VAE_func_DK(x0,
                        kn_method$latent_dim,
                        kn_method$lambda_kl,
                        kn_method$lambda_abun,
                        kn_method$lambda_pres,
                        kn_method$gamma_full,
                        kn_method$gamma_swap,
                        kn_method$lambda_moments,
                        kn_method$delta_corr,
                        kn_method$sigma_list,
                        kn_method$epochs,
                        kn_method$batch_size ,
                        kn_method$lr,
                        kn_method$weight_decay,
                        kn_method$seed)$knockoff_x
    }
    
  }else{
    
    x0 <- x
    
    if (kn_method$model == "2order"){
      xk <- knockoff::create.second_order(X=x0)
    }else if (kn_method$model == "VAE"){
      xk <- VAE_func(x0)$knockoff_x
    }else if (kn_method$model == "DK"){
      xk <- VAE_func_DK(x0,
                        kn_method$latent_dim,
                        kn_method$lambda_kl,
                        kn_method$lambda_abun,
                        kn_method$lambda_pres,
                        kn_method$gamma_full,
                        kn_method$gamma_swap,
                        kn_method$lambda_moments,
                        kn_method$delta_corr,
                        kn_method$sigma_list,
                        kn_method$epochs,
                        kn_method$batch_size ,
                        kn_method$lr,
                        kn_method$weight_decay,
                        kn_method$seed)$knockoff_x
    }
    
  }
  
  p <- ncol(cbind(x,xk))
  p0 <- ncol(x0)
  
  if (a > 0){
    lambda0 <- max(abs(t(sfun) %*% cbind(x,xk)))/(a*n)
  }else if (a == 0){
    lambda0 <- max(abs(t(sfun) %*% cbind(x,xk)))/(1e-3*n)
  }
  
  if (is.null(lambda.min.ratio)) lambda.min.ratio = ifelse(n < p, 1e-02, 1e-02)
  
  lambda <- 10^(seq(log10(lambda0),log10(lambda0*lambda.min.ratio),length.out=length.lambda))
  
  if (progress) cat("Algorithm running for full dataset: \n")
  
  fullfit <- logistic_enet_al_kn(cbind(x,xk),y,length.lambda,mu,maxiter,lambda,wcov,a,adjust,ncov,progress,FALSE,FALSE)
  
  if (!is.null(colnames(x))){
    rownames(fullfit$beta) = c(colnames(x),paste0(colnames(x)[(ncov+1):ncol(x)],"_kn"))
  }else{
    colnames(x) = c(paste0("cov",1:ncov),paste0("feat",1:p0))
    rownames(fullfit$beta) = c(colnames(x),paste0(colnames(x)[(ncov+1):ncol(x)],"_kn"))
  }
  
  df_betaseq <- data.frame(t(fullfit$beta[(ncov+1):p,]))
  df_betaseq$lambda <- lambda
  df_betaseq_long <- df_betaseq |> 
    pivot_longer(all_of(rownames(fullfit$beta[(ncov+1):p,])),
                 names_to = "feature",
                 values_to = "coef") |> 
    filter(coef != 0) |> 
    group_by(feature) |> 
    filter(row_number() == 1) |> 
    ungroup() |> 
    rowwise() |> 
    mutate(feature_original = strsplit(feature,"_")[[1]][1],
           kn = grepl("kn",feature)) |> 
    select(-coef,-feature) |> 
    pivot_wider(names_from = "kn",
                values_from = "lambda") |> 
    rename(KN = `TRUE`,
           Original = `FALSE`) 
  
  df_betaseq_long[is.na(df_betaseq_long)] <- 1e-8
  
  df_betaseq_long <- df_betaseq_long |> 
    mutate(W = ifelse(KN > Original,
                      -KN,
                      ifelse(KN < Original,
                             Original,0)))
  
  threshold <- knockoff::knockoff.threshold(df_betaseq_long$W,fdr=fdr,offset=kn_method$offset)
  thres.idx <- which(lambda == threshold)
  thres.beta = fullfit$beta[,thres.idx][(ncov+1):(ncov+p0)]
  selected.features = names(thres.beta[thres.beta != 0])
  
  ret <- list(intercept=fullfit$beta0,
              beta=fullfit$beta,
              lambda=fullfit$lambda,
              a=a,
              loss=fullfit$loss,
              mse=fullfit$mse,
              threshold=threshold,
              thres.beta = thres.beta,
              selected.features = selected.features
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
    
    pkn <- df_betaseq_long |> 
      ggplot(aes(x=Original,y=KN)) + 
      geom_text(aes(label = feature_original),alpha=1) +
      xlab("Lambda when original features entered") +
      ylab("Lambda when knockoff features entered") +
      geom_abline(slope=1,intercept = 0,color="red",linetype=2)+ 
      geom_vline(xintercept=threshold,color="red",linetype=2)
    
    ret$pkn <- pkn
    
    knockoffgp <- factor(c(rep("original",p0),rep("knockoff",p0)),
                         levels = c("original","knockoff"))
    
    pheatx <- suppressMessages(ComplexHeatmap::Heatmap (cbind(x0,xk),
                                                        name = "Log count",
                                                        cluster_columns = F,
                                                        cluster_rows = F,
                                                        show_row_names = F,
                                                        show_column_names = F,
                                                        column_split = knockoffgp))
    
    ret$data_matrix_heatmap <- pheatx
    
    pheatcov <- suppressMessages(ComplexHeatmap::Heatmap (var(cbind(x0,xk)),
                                                          name = "Covariance",
                                                          cluster_columns = F,
                                                          cluster_rows = F,
                                                          show_row_names = F,
                                                          show_column_names = F,
                                                          row_split = knockoffgp,
                                                          column_split = knockoffgp))
    
    ret$covariance_heatmap <- pheatcov
  }
  
  ret$runtime = proc.time() - ptm
  
  return(ret)
  
}