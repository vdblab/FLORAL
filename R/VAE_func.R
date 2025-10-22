VAE_func <- function(x,
                     latent_dim = 32,
                     lambda_kl = 1,
                     lambda_pres = 1){
  
  require(torch)
  require(luz)
  
  HybridVAE <- nn_module(
    "HybridVAE",
    initialize = function(input_dim, 
                          latent_dim, 
                          lambda_kl, 
                          lambda_pres) {
      self$encoder <- nn_sequential(
        nn_linear(input_dim, 256),
        nn_relu(),
        nn_linear(256, 128),
        nn_relu(),
        nn_linear(128, 64),
        nn_relu()
      )
      self$z_mean <- nn_linear(64, latent_dim)
      self$z_log_var <- nn_linear(64, latent_dim)
      self$lambda_kl = lambda_kl
      self$lambda_pres = lambda_pres
      
      self$decoder <- nn_sequential(
        nn_linear(latent_dim, 64),
        nn_relu(),
        nn_linear(64, 128),
        nn_relu(),
        nn_linear(128, 256),
        nn_relu(),
        nn_linear(256, input_dim)
      )
      
      self$presence_decoder <- nn_sequential(
        nn_linear(latent_dim, 64),
        nn_relu(),
        nn_linear(64, input_dim),
        nn_sigmoid()
      )
      
    },
    
    encode = function(x) {
      h <- self$encoder(x)
      mu <- self$z_mean(h)
      logvar <- self$z_log_var(h)
      list(mu, logvar)
    },
    
    reparameterize = function(mu, logvar) {
      std <- torch_exp(0.5 * logvar)
      eps <- torch_randn_like(std)
      mu + eps * std
    },
    
    decode = function(z) {
      self$decoder(z)
    },
    
    decode_presence = function(z) {
      self$presence_decoder(z)
    },
    
    forward = function(x) {
      encoded <- self$encode(x)
      mu <- encoded[[1]]
      logvar <- encoded[[2]]
      z <- self$reparameterize(mu, logvar)
      recon_x <- self$decode(z)
      pres_prob <- self$decode_presence(z)
      list(recon_x=recon_x, 
           input_x=x,
           pres_prob=pres_prob, 
           mu=mu, 
           logvar=logvar)
    },
    
    ### Loss function without incorporating deep knockoff
    
    loss = function(input, target){
      
      recon_x = input$recon_x
      input_x = torch_tensor(input$input_x, dtype = torch_float())
      pres_prob = input$pres_prob
      mu = input$mu
      logvar = input$logvar
      
      mask <- torch_tensor(input$input_x > 0, dtype = torch_float()) #(input_x > 0)$to(dtype = torch_float())
      
      bernoulli <- distr_bernoulli(probs = pres_prob)
      log_prob_presence <- bernoulli$log_prob(mask)
      loss_presence <- -log_prob_presence$sum()
      
      loss_abundance <- nnf_mse_loss(recon_x * mask, input_x, reduction = "sum")
      
      kl_loss <- -0.5 * torch_sum(1 + logvar - mu$pow(2) - torch_exp(logvar))
      
      batch_size <- input_x$size(1)
      loss <- (self$lambda_pres * loss_presence + loss_abundance + self$lambda_kl * kl_loss) / batch_size
      
      return(loss)
      
      # list(loss = loss,
      #      loss_presence = loss_presence,
      #      loss_abundance = loss_abundance,
      #      kl_loss = kl_loss)
      
    }
    
  )
  
  ### Setup a luz module to train the VAE model
  
  model <- HybridVAE %>%
    setup(
      loss = NULL,
      optimizer = optim_adamw
    ) %>%
    set_opt_hparams(weight_decay = 0.01) |>
    set_hparams(input_dim = ncol(x),
                latent_dim = latent_dim,
                lambda_kl = lambda_kl,
                lambda_pres = lambda_pres
    )
  
  ### Fit the model via luz
  
  fitted <- model |> fit(
    x,
    epochs = 100, # Subject to tuning?
    dataloader_options = list(batch_size = 50,shuffle=TRUE),
    accelerator = accelerator(cpu=TRUE),
    verbose = FALSE
  )
  
  ### Obtain the output from the trained model
  
  fitted$model$eval()
  
  with_no_grad({
    output <- fitted$model(torch_tensor(x,dtype=torch_float()))
  })
  
  input_x <- as_array(output$input_x)
  recon_x <- as_array(output$recon_x)
  pres_prob <- as_array(output$pres_prob)
  mask <- as_array(torch_bernoulli(output$pres_prob))
  
  ### Compute knockoff matrix
  knockoff_x <- recon_x * mask #(input_x > 0)
  
  return(list(knockoff_x = knockoff_x,
              recon_x = recon_x)
         )
  
}