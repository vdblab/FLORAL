simulation <- function(n = 100, p = 1000, true_size = 10){
  
  true_set <- 1:true_size
  beta <- rnorm(true_size,mean=0,sd=1)
  x <- matrix(NA,nrow=n,ncol=p)
  y <- rep(NA,length=n)
  for (j in 1:p){
    
    x[,j] <- scale(rnorm(n,mean=0,sd=1))
    # x[,j] <- scale(rmvnorm(1,rep(0,p),sigma=diag(p))) 
    
  }
  
  y <- x[,true_set] %*% beta + rnorm(n,mean=0,sd=1)
  y0 = y - mean(y)
  
  return(list(x=x,y=y,y0=y0,beta=beta,idx=true_set))
  
}

simulation_constrained <- function(n = 100, p = 1000, true_size = 10, intercept=FALSE){
  
  true_set <- 1:true_size
  beta <- c(1,-1,1,2,-1,1,-1,-1,1)
  beta[true_size] = 0-sum(beta[-true_size])
  x <- matrix(NA,nrow=n,ncol=p)
  y <- rep(NA,length=n)
  for (j in 1:p){
    
    x[,j] <- scale(rnorm(n,mean=0,sd=1))
    # x[,j] <- scale(rmvnorm(1,rep(0,p),sigma=diag(p))) 
    
  }
  
  y <- x[,true_set] %*% beta + rnorm(n,mean=0,sd=0.1)
  
  if(intercept) {
    
    intcpt <- rnorm(1,mean=1,sd=1)
    y <- y + intcpt
    
  }
  y0 = y - mean(y)
  
  ret <- list(x=x,y=y,y0=y0,beta=beta,idx=true_set)
  
  if (intercept) ret$intercept=intcpt
  
  return(ret)
  
}

simulation_constrained_binary <- function(n = 100, p = 1000, true_size = 10){
  
  true_set <- 1:true_size
  beta <- c(2,2,2,2,2,-2,-2,-2,-2)
  beta[true_size] = 0-sum(beta[-true_size])
  x <- matrix(NA,nrow=n,ncol=p)
  y <- rep(NA,length=n)
  for (j in 1:p){
    
    x[,j] <- scale(rnorm(n,mean=0,sd=1))
    # x[,j] <- scale(rmvnorm(1,rep(0,p),sigma=diag(p))) 
    
  }
  
  eta <- x[,true_set] %*% beta
  
  prob <- exp(eta)/(1+exp(eta))
  
  for (i in 1:n){
    y[i] <- rbinom(1,1,prob=prob[i])
  }
  
  return(list(x=x,y=y,beta=beta,idx=true_set))
  
}


simulation_constrained_cox <- function(n = 100, p = 1000, true_size = 10){
  
  true_set <- 1:true_size
  beta <- c(2,2,2,2,2,-2,-2,-2,-2)
  beta[true_size] = 0-sum(beta[-true_size])
  x <- matrix(NA,nrow=n,ncol=p)
  t0 <- rep(NA,length=n)
  for (j in 1:p){
    
    x[,j] <- scale(rnorm(n,mean=0,sd=1))
    # x[,j] <- scale(rmvnorm(1,rep(0,p),sigma=diag(p))) 
    
  }
  
  eta <- x[,true_set] %*% beta
  
  lambda <- exp(eta)
  
  for (i in 1:n){
    t0[i] <- rexp(1,rate=lambda[i])
  }
  
  c <- runif(n,min=0,max=10)
  
  t <- ifelse(t0 <= c, t0, c)
  d <- as.numeric(t0 <= c)
  
  return(list(x=x,t=t,d=d,beta=beta,idx=true_set))
  
}
