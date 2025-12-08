LogRatioGEEKN <- function(x,
                          y,
                          id,
                          ncov,
                          intercept,
                          family,
                          corstr,
                          scalefix,
                          scalevalue,
                          length.lambda=100,
                          lambda.min.ratio=NULL,
                          wcov,
                          a=1,
                          mu=1e6,
                          pfilter=0,
                          fdr=0.05,
                          kn_method=NULL,
                          offset=c(0,1),
                          maxiter1=100,
                          maxiter2=1,
                          progress=TRUE,
                          plot=TRUE){
  
  ptm <- proc.time()
  
  n <- length(unique(id))
  p <- ncol(x)
  
  adjust = FALSE
  if (ncov > 0) adjust = TRUE
  
  if (adjust){
    
    x0 <- x[,(ncov+1):ncol(x)]
    
    if (kn_method$model == "2order"){
      xk <- knockoff::create.second_order(X=x0)
    }else if (kn_method$model == "VAE"){
      xk <- train_vae(x0,
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
                      kn_method$seed,
                      progress=progress)$knockoff_x
    }
    
  }else{
    
    x0 <- x
    
    if (kn_method$model == "2order"){
      xk <- knockoff::create.second_order(X=x0)
    }else if (kn_method$model == "VAE"){
      xk <- train_vae(x0,
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
                      kn_method$seed,
                      progress=progress)$knockoff_x
    }
    
  }
  
  p_total <- ncol(cbind(x,xk))
  p0 <- ncol(x0)
  
  if (family == "gaussian"){
    
    if (a > 0){
      lambda0 <- max(abs(t(scale(y)) %*% cbind(x,xk)))/(min(a,1)*nrow(cbind(x,xk)))*100
    }else if (a == 0){
      lambda0 <- max(abs(t(scale(y)) %*% cbind(x,xk)))/(1e-3*nrow(cbind(x,xk)))*100
    }
    
  }else if (family == "binomial"){
    
    sfun = y-0.5
    
    if (a > 0){
      lambda0 <- max(abs(t(sfun) %*% cbind(x,xk)))/(min(a,1)*nrow(cbind(x,xk)))*100
    }else if (a == 0){
      lambda0 <- max(abs(t(sfun) %*% cbind(x,xk)))/(1e-3*nrow(cbind(x,xk)))*100
    }
    
  }
  
  family0 <- family
  if (is.character(family)) family <- get(family)
  if (is.function(family))  family <- family()
  
  if (is.null(lambda.min.ratio)) lambda.min.ratio = ifelse(n < p_total, 1e-01, 1e-02)
  lambda <- 10^(seq(log10(lambda0),log10(lambda0*lambda.min.ratio),length.out=length.lambda))
  
  nt <- as.integer(unlist(lapply(split(id, id), "length"))) # Number of obs per subject
  
  if (progress) cat("Algorithm running for full dataset: \n")
  
  fullfit <- gee_fit_kn(y,
                        cbind(x,xk),
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
                        maxiter1=maxiter1,
                        maxiter2=maxiter2,
                        scalefix=scalefix,
                        scalevalue=scalevalue,
                        display_progress=progress)
  
  if (!is.null(colnames(x))){
    rownames(fullfit$beta) = c(colnames(x),paste0(colnames(x)[(ncov+1):ncol(x)],"_kn"))
  }else{
    colnames(x) = c(paste0("cov",1:ncov),paste0("feat",1:p0))
    rownames(fullfit$beta) = c(colnames(x),paste0(colnames(x)[(ncov+1):ncol(x)],"_kn"))
  }
  
  beta_filtered <- fullfit$beta
  beta_filtered[abs(beta_filtered) < 1e-3] = 0
  
  if (ncov > 0){
    idxfeat <- setdiff(1:nrow(beta_filtered),1:ncov)
  }else{
    idxfeat <- 1:nrow(beta_filtered)
  }
  
  beta_filtered[idxfeat,][apply(beta_filtered[idxfeat,], 2,function(x) abs(x) < max(abs(x))*pfilter)] <- 0
  
  df_betaseq <- data.frame(t(beta_filtered[(ncov+1):p_total,]),check.names = FALSE)
  df_betaseq$lambda <- lambda
  df_betaseq_long <- df_betaseq |> 
    pivot_longer(all_of(rownames(beta_filtered[(ncov+1):p_total,])),
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
  
  # Ensure fdr and offset are vectors
  fdr <- as.vector(fdr)
  offset <- as.vector(offset)
  
  # Calculate threshold, thres.idx, thres.beta, and selected.features for each FDR×offset combination
  threshold_list <- list()
  thres.idx_list <- list()
  thres.beta_list <- list()
  selected.features_list <- list()
  
  combo_idx <- 1
  for (i in seq_along(fdr)) {
    for (j in seq_along(offset)) {
      fdr_val <- fdr[i]
      offset_val <- offset[j]
      threshold <- knockoff::knockoff.threshold(df_betaseq_long$W, fdr=fdr_val, offset=offset_val)
      threshold_list[[combo_idx]] <- threshold
      names(threshold_list)[combo_idx] <- paste0("fdr_", fdr_val, "_offset_", offset_val)
      
      # Find the lambda index closest to the threshold (in case of floating point precision issues)
      thres.idx <- which.min(abs(lambda - threshold))
      thres.idx_list[[combo_idx]] <- thres.idx
      names(thres.idx_list)[combo_idx] <- paste0("fdr_", fdr_val, "_offset_", offset_val)
      
      thres.beta <- beta_filtered[,thres.idx][(ncov+1):(ncov+p0)]
      thres.beta_list[[combo_idx]] <- thres.beta
      names(thres.beta_list)[combo_idx] <- paste0("fdr_", fdr_val, "_offset_", offset_val)
      
      selected.features <- names(thres.beta[thres.beta != 0])
      selected.features_list[[combo_idx]] <- selected.features
      names(selected.features_list)[combo_idx] <- paste0("fdr_", fdr_val, "_offset_", offset_val)
      
      combo_idx <- combo_idx + 1
    }
  }
  
  # For backward compatibility, also store the first combination as scalar outputs
  threshold <- threshold_list[[1]]
  thres.idx <- thres.idx_list[[1]]
  thres.beta <- thres.beta_list[[1]]
  selected.features <- selected.features_list[[1]]
  
  ret <- list(beta=beta_filtered,
              beta_unfiltered = fullfit$beta,
              tol=fullfit$tol,
              iters=fullfit$iters,
              lambda=lambda,
              a=a,
              threshold=threshold,  # First combination for backward compatibility
              thres.beta = thres.beta,  # First combination for backward compatibility
              selected.features = selected.features,  # First combination for backward compatibility
              fdr = fdr,  # Store the FDR vector
              offset = offset,  # Store the offset vector
              threshold_list = threshold_list,  # All thresholds (all FDR×offset combinations)
              thres.idx_list = thres.idx_list,  # All threshold indices
              thres.beta_list = thres.beta_list,  # All threshold betas
              selected.features_list = selected.features_list  # All selected features
  )
  
  if (plot){
    
    beta_nzero <- suppressWarnings(data.frame(reshape::melt(ret$beta[rowSums(ret$beta != 0) > 0,])))
    beta_nzero$lambda <- ret$lambda[beta_nzero$X2]
    beta_nzero$loglambda <- log(beta_nzero$lambda)
    
    lambda_count <- data.frame(loglambda = log(ret$lambda),
                               count = colSums(ret$beta != 0))
    lambda_count <- lambda_count[seq(5,nrow(lambda_count),length.out=10),]
    
    top10feat <- sort(ret$beta[,length(ret$lambda)])[c(1:5,(p_total-4):p_total)]
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
    
    # Create base plot
    pkn <- df_betaseq_long |> 
      ggplot(aes(x=Original,y=KN)) + 
      geom_text(aes(label = feature_original),alpha=1) +
      xlab("Lambda when original features entered") +
      ylab("Lambda when knockoff features entered") +
      geom_abline(slope=1,intercept = 0,color="red",linetype=2)
    
    # Add vertical lines for each FDR×offset threshold
    if (length(threshold_list) > 1) {
      # Multiple combinations - use different colors
      colors <- rainbow(length(threshold_list))
      for (i in seq_along(threshold_list)) {
        threshold_val <- threshold_list[[i]]
        combo_name <- names(threshold_list)[i]
        # Extract fdr and offset from name (format: "fdr_X_offset_Y")
        parts <- strsplit(combo_name, "_")[[1]]
        fdr_val <- parts[2]
        offset_val <- parts[4]
        pkn <- pkn + 
          geom_vline(xintercept=threshold_val, 
                     color=colors[i], 
                     linetype=2, 
                     alpha=0.8,
                     linewidth=1) +
          annotate("text", 
                   x=threshold_val, 
                   y=max(df_betaseq_long$KN, na.rm=TRUE) * (0.95 - (i-1)*0.05), 
                   label=paste0("FDR=", fdr_val, "\nOff=", offset_val), 
                   color=colors[i],
                   angle=90,
                   vjust=-0.5,
                   hjust=0.5,
                   size=2.5)
      }
    } else {
      # Single threshold case - use red
      pkn <- pkn + 
        geom_vline(xintercept=threshold, color="red", linetype=2)
    }
    
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
  
  ret$runtime <- proc.time() - ptm
  
  return(ret)
}

