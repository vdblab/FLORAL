gee_NR_R <-
  function(N, # Number of subjects
           nt, # number of obs per subject
           y,
           X,
           nx, # Number of covariates (ncol(X))
           family,
           beta_new,
           Rhat, # estimated working correlation matrix (Ehat from cor_gee.R)
           fihat, # estimated scale parameter (fi from cor_gee.R)
           lambda,
           pindex, # index of parameters not subject to penalization (!!!!!! Can create another index for zero-sum constraint)
           eps=10^-6,
           muu) {
    
    aindex=cumsum(nt)
    index=c(0,aindex[-length(aindex)])
    
    eta=X%*%beta_new
    mu=family$linkinv(eta)
    
    #This is E on Wang et al.(2012)
    E1 <- diag(q_lasso(abs(as.vector(beta_new)),lambda)/(abs(as.vector(beta_new))+eps)) 
    
    # E2 <- diag(x=muu*(sum(as.vector(beta_new))-alpha),nrow=length(as.vector(beta_new)),ncol=length(as.vector(beta_new)))
    E2 <- rep(muu*(sum(as.vector(beta_new))),length(as.vector(beta_new)))
    # print(dim(E1))
    # print(dim(E2))
    
    if(is.null(pindex)==TRUE) {
      E<-E1 
    } else 
      if(is.null(pindex)!=TRUE) {
        E1[,pindex]<-0
        E<-E1
      }
    
    sum201<-matrix(0,nx,1)      #gradient:S
    sum301<-matrix(0,nx,nx)     #naive variance:H
    sum401<-matrix(0,nx,nx)     #a component for robust variance:M
    
    for (i in 1:N) {
      ym<-matrix(0,nt[i],1)
      bigD<-matrix(0,nt[i],nx)
      bigA<-matrix(0,nt[i],nt[i])
      for (j in 1:nt[i]) {
        #cat("j",j,"\n")
        ym[j]<- y[j+index[i]]-mu[j+index[i]] 
        bigA[j,j]<-family$variance(mu)[j+index[i]]
        for (k in 1:nx) {
          bigD[j,k]<-family$mu.eta(eta)[j+index[i]]*X[j+index[i],k]
          #cat("i",i,"j",j,"k",k,"\n")
        } # for k
      } # for j
      
      ##working covariance matrix
      bigV<-sqrt(bigA)%*%Rhat[1:nt[i],1:nt[i],i]%*%sqrt(bigA)
      #bigV<-fihat*bigV
      
      ##This is S in Wang et al.(2012)
      sum200<-t(bigD)%*%ginv(bigV)%*%ym      #this is like gradient
      sum201<-sum201+sum200
      
      ##This is H in Wang et al.(2012)
      sum300<-t(bigD)%*%ginv(bigV)%*%bigD    #this is for information matrix.
      sum301<-sum301+sum300
      
      ##Speed up the code##
      SSA=sqrt(ginv(bigA))
      SRhat=ginv(Rhat[1:nt[i],1:nt[i],i])
      SSAym=(SSA%*%ym)
      
      sum400<-t(bigD)%*%SSA%*%SRhat%*%(SSAym%*%t(SSAym))%*%SRhat%*%SSA%*%bigD
      sum401<-sum401+sum400
      
      #cat("i",i,"sum201",sum201,"sum301",sum301,"sum401",sum401,"\n")
    }  #end of i
    
    S<-fihat*sum201
    H<-fihat*sum301
    E<-E
    M<-fihat*sum401
    
    return(list("S"=S,"H"=H,"E"=E,"M"=M,"E2"=E2))
  }

q_scad <-
  function(theta,lambda,a=3.7)
  {
    #length of parameter
    p<-length(theta)
    #get the absolute value
    theta<-abs(theta)
    #create vector of zeros
    b1<-rep(0,p)
    #if theta is greater then lambda set it to 1
    b1[theta>lambda]<-1
    #create an another vector of zeros
    b2<-rep(0,p)
    #if theta is less than a*lambda, set it to 1.
    b2[theta<(lambda*a)]<-1
    lambda*(1-b1)+((lambda*a)-theta)*b2/(a-1)*b1
  }

q_lasso <- 
  function(theta,lambda){
    p <- length(theta)
    rep(lambda,p)
  }
