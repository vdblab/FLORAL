VAE_func_DK <- function(x,
                        latent_dim      = 32,
                        lambda_kl       = 0,
                        lambda_abun = 0,
                        lambda_pres     = 0,
                        gamma_full       = 1,
                        gamma_swap       = 1,
                        lambda_moments  = 0,
                        delta_corr      = 0,
                        sigma_list      = c(1,2,4,8,16,32,64,128),
                        epochs          = 100,
                        batch_size      = 50,
                        lr              = 1e-3,
                        weight_decay    = 1e-2,
                        seed            = 123) {
  
  require(torch); require(luz)
  torch_manual_seed(seed); set.seed(seed)
  
  p          <- ncol(x)
  SigmaHat_t <- torch_tensor(cov(x),              dtype = torch_float()) # covariance matrix of x
  Mask_t     <- torch_tensor(matrix(1, p, p)-diag(p), dtype = torch_float()) # mask matrix: keep off-diagonal entries only
  target_t   <- torch_zeros(p)
  
  mix_rbf_mmd2_loss <- function(X, Y, sigmas) {
    XX <- torch_cdist(X, X, p = 2)$pow(2) # 2-norm distances ^ 2
    YY <- torch_cdist(Y, Y, p = 2)$pow(2)
    XY <- torch_cdist(X, Y, p = 2)$pow(2)
    
    mmd2 <- torch_tensor(0, dtype = torch_float(), device = X$device)
    
    for (s in sigmas) {
      g <- 1 / (2*s^2) # s is the bandwidth parameter
      mmd2 <- mmd2 +
        torch_mean(torch_exp(-g*XX)) + # Gaussian kernel
        torch_mean(torch_exp(-g*YY)) -
        2*torch_mean(torch_exp(-g*XY))
    }
    torch_sqrt(torch_relu(mmd2))
  }
  
  HybridVAE <- nn_module(
    "HybridVAE",
    initialize = function(input_dim) {
      self$enc <- nn_sequential(
        nn_linear(input_dim,256),nn_relu(),
        nn_linear(256,128),nn_relu(),
        nn_linear(128,64), nn_relu())
      self$z_mu     <- nn_linear(64, latent_dim)
      self$z_logvar <- nn_linear(64, latent_dim)
      
      self$dec <- nn_sequential(
        nn_linear(latent_dim,64),nn_relu(),
        nn_linear(64,128),nn_relu(),
        nn_linear(128,256),nn_relu(),
        nn_linear(256,input_dim))
      self$pres_dec <- nn_sequential(
        nn_linear(latent_dim,64),nn_relu(),
        nn_linear(64,input_dim),nn_sigmoid())
    },
    
    encode = function(x){
      h <- self$enc(x)
      list(self$z_mu(h), self$z_logvar(h))
    },
    reparam = function(mu, lv) mu + torch_randn_like(mu)*torch_exp(0.5*lv),
    decode      = function(z) self$dec(z),
    decode_pres = function(z) self$pres_dec(z),
    
    forward = function(x){
      mu_lv <- self$encode(x)
      z     <- self$reparam(mu_lv[[1]], mu_lv[[2]])
      list(
        recon_x   = self$decode(z),
        pres_prob = self$decode_pres(z),
        mu        = mu_lv[[1]],
        logvar    = mu_lv[[2]],
        input_x   = x)
    },
    
    loss = function(input, target){
      
      recon_x   <- input$recon_x
      input_x   <- torch_tensor(input$input_x, dtype = torch_float())
      pres_prob <- input$pres_prob
      mu        <- input$mu
      logvar    <- input$logvar
      
      mask      <- torch_tensor(input$input_x > 0, dtype = torch_float())
      n_batch   <- input_x$size(1)
      
      bern      <- distr_bernoulli(probs = pres_prob)
      loss_pres <- -bern$log_prob(mask)$sum()
      
      xk <- recon_x * pres_prob
      loss_abun <- nnf_mse_loss(xk, input_x, reduction = "sum")
      kl_loss   <- -0.5 * torch_sum(1 + logvar - mu$pow(2) - torch_exp(logvar))
      
      n_half <- n_batch %/% 2
      X1  <- input_x[1:n_half, ];     X2  <- input_x[(n_half+1):(2*n_half), ]
      Xk1 <- xk[1:n_half, ];     Xk2 <- xk[(n_half+1):(2*n_half), ]
      
      Z1 <- torch_cat(list(X1,Xk1), dim = 2)
      Z2 <- torch_cat(list(Xk2,X2), dim = 2) # Changed here: wrong definition
      Z3 <- torch_cat(list(X2,Xk2), dim = 2)
      p  <- input_x$size(2)
      swap <- which(rbinom(p, 1, 0.5) == 1)
      if (length(swap)){
        Z3[, swap]     <- Xk2[, swap]
        Z3[, swap + p] <- X2[, swap]
      }
      mmd_full <- mix_rbf_mmd2_loss(Z1, Z2, sigma_list)
      mmd_swap <- mix_rbf_mmd2_loss(Z1, Z3, sigma_list)
      
      mX  <- input_x - input_x$mean(1, keepdim = TRUE)
      mXk <- xk - xk$mean(1, keepdim = TRUE)
      SXk  <- torch_matmul(mXk$t(), mXk)/n_batch
      SXXk <- torch_matmul(mX$t(),  mXk)/n_batch
      loss_1m <- (input_x$mean(1) - xk$mean(1))$pow(2)$sum()
      scale_f <- torch_sum(SigmaHat_t$pow(2))
      loss_2m <- ((SigmaHat_t - SXk)$pow(2)$sum() +
                    (Mask_t * (SigmaHat_t - SXXk))$pow(2)$sum()) / scale_f
      loss_mom <- loss_1m + loss_2m
      
      scaleX  <- mX$pow(2)$mean(1, keepdim = TRUE)
      scaleXk <- mXk$pow(2)$mean(1, keepdim = TRUE)
      corr <- ((mX  / (torch_sqrt(scaleX)  + 1e-3)) * # location of eps changed
                 (mXk / (torch_sqrt(scaleXk) + 1e-3)))$mean(1)
      loss_corr <- (corr - target_t)$pow(2)$mean()
      
      loss <- (lambda_pres    * loss_pres +
          lambda_abun * loss_abun +
          lambda_kl      * kl_loss +
          gamma_full      * mmd_full + 
          gamma_swap      * mmd_swap +
          lambda_moments * loss_mom +
          delta_corr     * loss_corr) / n_batch
      
      return(loss)
    }
  )
  
  model <- HybridVAE %>% setup(loss = NULL, optimizer = optim_adamw) %>%
    set_hparams(input_dim = ncol(x)) %>%
    set_opt_hparams(lr = lr, weight_decay = weight_decay)
  
  fitted <- model %>% fit(
    x,
    epochs = epochs,
    dataloader_options = list(batch_size = batch_size, shuffle = TRUE),
    accelerator = accelerator(cpu = TRUE))
  
  fitted$model$eval()
  with_no_grad({
    output <- fitted$model(torch_tensor(x, dtype = torch_float()))
  })
  
  input_x <- as_array(output$input_x)
  recon_x <- as_array(output$recon_x)
  pres_prob <- as_array(output$pres_prob)
  mask <- as_array(torch_bernoulli(output$pres_prob))
  
  knockoff_x <- recon_x * mask #(input_x > 0)
  
  return(list(knockoff_x = knockoff_x,
              recon_x = recon_x)
  )
}
