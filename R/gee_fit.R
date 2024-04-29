gee_fit <- function(x,
                    y,
                    id,
                    family,
                    corstr,
                    lambda,
                    maxiter=100,
                    tol=1e-4,
                    scale.fix=0,
                    scale.value=NULL,
                    R=NULL,
                    pindex=NULL, # index of parameters not subject to penalization (!!!!!! Can create another index for zero-sum constraint)
                    eps=1e-6,
                    muu=1e6){
  
  N <- length(unique(id))
  nt <- as.integer(unlist(lapply(split(id, id), "length"))) # Number of obs per subject
  maxclsz <-max(nt) # Max number of obs
  nobs <- length(y)
  p <- ncol(x)
  
  beta_mat <- matrix(NA,nrow=p,ncol=length(lambda))
  
  beta_new = rep(0,p)
  
  for (j in 1:length(lambda)){
    
    for (i in 1:maxiter){
      
      beta <- beta_new
      
      # cor.obj <- gee_cor_R(N,
      #                      nt,
      #                      y,
      #                      x,
      #                      family, 
      #                      beta,
      #                      corstr,
      #                      Mv=NULL,
      #                      maxclsz, # max number of obs
      #                      R, # User-specified structure
      #                      scale.fix, # 0 or 1: indicator of fixed scale parameter
      #                      scale.value # Value of the scale parameter (if fixed)
      # )
      # Rhat=cor.obj$Rhat
      # fihat=cor.obj$fi
      
      
      cor.obj = gee_cor(
        N,
        nt,
        y,
        x,
        family$linkinv,
        family$variance,
        beta,
        corstr,
        maxclsz, # max number of obs
        scale.fix==0,
        1
      )
      
      Rhat=cor.obj$Rhat
      fihat=cor.obj$fihat
      
      # NR.obj=gee_NR_R(N,
      #               nt,
      #               y,
      #               x,
      #               p,
      #               family,
      #               beta,
      #               Rhat,
      #               fihat,
      #               lambda[j],
      #               pindex,
      #               eps,
      #               muu)
      # S<-NR.obj$S
      # H<-NR.obj$H
      # E<-NR.obj$E
      # C<- NR.obj$E2
      
      NR.obj=gee_NR(N,
                    nt,
                    y,
                    x,
                    p,
                    family$linkinv,
                    family$mu.eta,
                    family$variance,
                    beta_new,
                    Rhat,
                    fihat,
                    lambda[j],
                    eps,
                    muu)
      
      S<-NR.obj$S
      H<-NR.obj$H
      E<-NR.obj$E
      C<- NR.obj$C
      
      beta_new <- matrix(beta)+ginv(H+N*E+muu)%*%(S-N*E%*%matrix(beta)-C)
      
      diff<-sum(abs(beta-beta_new)) 
      
      if (verbose) print(paste0("Lambda:",round(lambda[j],3),", iter:",i,", diff:",round(diff,6)))
      if (diff <= tol) break
      
    }
    
    beta_mat[,j] <- beta_new
    
  }
  
  beta_mat_raw <- beta_mat
  beta_mat[abs(beta_mat) < 1e-3] <- 0
  
}