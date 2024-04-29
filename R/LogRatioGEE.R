LogRatioGEE <- function(x,
                        y,
                        id,
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
  
  if (a > 0){
    lambda0 <- max(abs(t(y) %*% x))/(a*n)
  }else if (a == 0){
    lambda0 <- max(abs(t(y) %*% x))/(1e-3*n)
  }
  
  if (is.null(lambda.min.ratio)) lambda.min.ratio = ifelse(n < p, 1e-02, 1e-02)
  lambda <- 10^(seq(log10(lambda0),log10(lambda0*lambda.min.ratio),length.out=length.lambda))
  
  if (progress) cat("Algorithm running for full dataset: \n")
  
  fullfit <- gee_fit(x,y,id,length.lambda,mu,100,lambda,wcov,a,adjust,ncov,progress)
  
  
  
  if(!(is.double(x)))  x <- as.double(x)
  if(!(is.double(y)))  y <- as.double(y)
  if(!(is.double(id))) id <- as.double(id)
  
  N<-length(unique(id))
  K=ncol(x)-1
  nx=ncol(x)
  
  avec <- as.integer(unlist(lapply(split(id, id), "length")))
  maxclsz <-max(avec)
  maxcl <- maxclsz
  nt<-avec
  nobs<-sum(nt)
  
  xnames <- dimnames(X)[[2]]
  if(is.null(xnames) && colnames(X)[1]=="(Intercept)") {
    xnames <- paste("x", 0:K, sep = "")
    dimnames(X) <- list(NULL, xnames)
  } else
    if(is.null(xnames) && colnames(X)[1]!="(Intercept)") {
      xnames <- paste("x", 1:K, sep = "")
      dimnames(X) <- list(NULL, xnames)
    }
  
  if(!(is.double(N)))      N <- as.double(N)
  if(!(is.double(maxcl)))  maxcl <- as.double(maxcl)
  if(!(is.double(nobs)))   nobs <- as.double(nobs)
  
  if(missing(lambda)) stop("A value is not assiged for lambda!")
  
  if(missing(muu)) stop("A value is not assiged for muu!")
  
  if(missing(pindex)) pindex=NULL
  
  if(missing(family)) family=gaussian(link="identity")
  
  if(missing(corstr)) corstr="independence"
  
  if(missing(Mv)) Mv<-NULL
  
  if(corstr=="stat_M_dep" && is.null(Mv)) stop("corstr is assumed to be 'stat_M_dep' but Mv is not specified!")
  
  if(corstr=="non_stat_M_dep" && is.null(Mv)) stop("corstr is assumed to be 'non_stat_M_dep' but Mv is not specified!")
  
  if((corstr!="stat_M_dep" && corstr!="non_stat_M_dep") && !is.null(Mv))  stop("Mv is specified while corstr is assumed to be neither 
'stat_M_dep' nor 'non_stat_M_dep'!")
  
  if(corstr=="non_stat_M_dep" && length(unique(nt)) !=1) stop("corstr cannot be assumed to be 'non_stat_M_dep' for unbalanced data!")
  
  if(corstr=="unstructured" && length(unique(nt)) !=1) stop("corstr cannot be assumed to be 'unstructured' for unbalanced data!")
  
  if(missing(R)) R<-NULL
  
  if(corstr=="fixed" && is.null(R))  stop("corstr is assumed to be 'fixed' but R is not specified!")
  if(corstr!="fixed" && !is.null(R)) stop("R is specified although corstr is not assumed to be 'fixed'!")
  
  if(!is.null(R)) {
    Rr <- nrow(R)
    if(Rr != ncol(R)) stop("R is not square!")
    if(Rr < maxclsz)  {stop("R is not big enough to accommodate some clusters!")} else
      if(Rr > maxclsz)  {stop("R is larger than the maximum cluster!")}
  }
  
  if(missing(scale.fix))  scale.fix <- TRUE
  scale.fix <- as.integer(scale.fix)
  
  if(missing(scale.value)) scale.value=1
  scale.value<-as.integer(scale.value)
  
  if(missing(eps)) eps=10^-6
  eps<-as.double(eps)
  
  if(missing(maxiter)) maxiter<-30
  maxiter<-as.integer(maxiter)
  
  if(missing(tol))  tol=10^-3
  tol=as.double(tol)
  
  if(missing(silent))  silent <-TRUE
  silent<-as.integer(silent)
  
  if (is.character(family)) family <- get(family)
  if (is.function(family))  family <- family()
  
  links <- c("identity","log","logit","inverse","probit","cloglog")
  fams <- c("gaussian","poisson","binomial","Gamma","quasi")
  varfuns <- c("constant", "mu", "mu(1-mu)", "mu^2")
  corstrs <- c("independence", "fixed", "stat_M_dep", "non_stat_M_dep", "exchangeable", 
               "AR-1", "unstructured")
  
  linkv <- as.integer(match(c(family$link), links, -1))
  if(linkv < 1) stop("unknown link!")
  
  famv <- match(family$family, fams, -1)
  if(famv < 1) stop("unknown family")
  if(famv <= 4) varfunv <- famv
  else varfunv <- match(family$varfun, varfuns, -1)
  if(varfunv < 1) stop("unknown varfun!")
  
  corstrv <- as.integer(match(corstr, corstrs, -1))
  if(corstrv < 1) stop("unknown corstr!")
  
  Mv <- as.integer(Mv)
  
  if (!is.null(beta_int))
  {
    beta <- matrix(beta_int, ncol = 1)
    if(nrow(beta) != nx) {stop("Dimension of beta != ncol(X)!")}
    #message("user\'s initial regression estimate")
    
  }
  else {
    #message("running glm to get initial regression estimate!")
    ### <tsl>	beta <- as.numeric(glm(m, family = family)$coef)
    mm <- match.call(expand.dots = FALSE)
    mm$R <- mm$beta_int <- mm$tol <- mm$maxiter <- mm$link <- 
      mm$varfun <-mm$corstr <- mm$Mv <- mm$silent <-mm$scale.fix <- 
      mm$scale.value <- mm$id<-
      mm$lambda <-mm$pindex<-mm$eps<-NULL
    mm[[1]]<-as.name("glm")
    beta <- eval(mm, parent.frame())$coef
    ### </tsl>
    print(beta)
    
  }
  
  # print(beta)
  
  beta_int=matrix(beta, ncol = 1)
  beta_new<-beta_int
  # muu <- 10
  alpha <- 0
  
  R.fi.hat=gee_cor(N,nt,y,X,family,beta_new,corstr,Mv,maxclsz,R=R,scale.fix=scale.fix,scale.value=scale.fix)
  Rhat=R.fi.hat$Ehat
  fihat=R.fi.hat$fi
  
  S.H.E.val=gee_NR(N,nt,y,X,nx,family,beta_new,Rhat,fihat,lambda,pindex,eps,muu=muu,alpha=alpha)
  S<-S.H.E.val$S
  H<-S.H.E.val$H
  E<-S.H.E.val$E
  E2 <- S.H.E.val$E2
  
  diff<-1
  iter<-0
  
  # while (iter_outer < maxiterouter){
  
  # beta_old_outer <- beta_new
  
  while(iter < maxiter) {
    
    beta_old<-beta_new
    
    beta_new<-matrix(beta_old)+ginv(H+N*E+muu)%*%(S-N*E%*%matrix(beta_old)-E2)
    
    R.fi.hat=gee_cor(N,nt,y,X,family,beta_new,corstr,Mv,maxclsz,R,scale.fix,scale.value)
    Rhat=R.fi.hat$Ehat
    fihat=R.fi.hat$fi
    
    # print(sum(beta_new))
    # print(alpha)
    
    S.H.E.M.val=gee_NR(N,nt,y,X,nx,family,beta_new,Rhat,fihat,lambda,pindex,eps,muu=muu,alpha=alpha)
    S<-S.H.E.M.val$S
    H<-S.H.E.M.val$H
    E<-S.H.E.M.val$E
    M<-S.H.E.M.val$M
    E2 <- S.H.E.M.val$E2
    
    diff<-sum(abs(beta_old-beta_new)) 
    
    iter<-iter+1
    # alpha <- alpha + sum(beta_new)
    if (silent==0) cat("iter",iter,"beta_new",beta_new,"diff",diff,"\n")
    if (diff <= tol) break
  } #end of while
  
  
}