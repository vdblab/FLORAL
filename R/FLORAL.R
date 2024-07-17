#' Fit Log-ratio lasso regression for compositional covariates 
#'
#' @description Conduct log-ratio lasso regression for continuous, binary and survival outcomes. 
#' @param x Feature matrix, where rows specify subjects and columns specify features. The first \code{ncov} columns should be patient characteristics and the rest columns are microbiome absolute counts corresponding to various taxa. If \code{x} contains longitudinal data, the rows must be sorted in the same order of the subject IDs used in \code{y}.
#' @param y Outcome. For a continuous or binary outcome, \code{y} is a vector. For survival outcome, \code{y} is a \code{Surv} object.
#' @param ncov An integer indicating the number of first \code{ncov} columns in \code{x} that will not be subject to the zero-sum constraint.
#' @param family Available options are \code{gaussian}, \code{binomial}, \code{cox}, \code{finegray}.
#' @param longitudinal \code{TRUE} or \code{FALSE}, indicating whether longitudinal data matrix is specified for input \code{x}. (\code{Longitudinal=TRUE} and \code{family="cox"} or \code{"finegray"} will fit a time-dependent covariate model. \code{Longitudinal=TRUE} and \code{family="gaussian"} or \code{"binomial"} will fit a GEE model.)
#' @param id If \code{longitudinal} is \code{TRUE}, \code{id} specifies subject IDs corresponding to the rows of input \code{x}.
#' @param tobs If \code{longitudinal} is \code{TRUE}, \code{tobs} specifies time points corresponding to the rows of input \code{x}.
#' @param failcode If \code{family = finegray}, \code{failcode} specifies the failure type of interest. This must be a positive integer.
#' @param corstr If a GEE model is specified, then \code{corstr} is the corresponding working correlation structure. Options are \code{independence}, \code{exchangeable}, \code{AR-1} and \code{unstructured}.
#' @param scalefix \code{TRUE} or \code{FALSE}, indicating whether the scale parameter is estimated or fixed if a GEE model is specified.
#' @param scalevalue Specify the scale parameter if \code{scalefix=TRUE}.
#' @param pseudo Pseudo count to be added to \code{x} before taking log-transformation. If unspecified, then the log-transformation will not be performed.
#' @param length.lambda Number of penalty parameters used in the path
#' @param lambda.min.ratio Ratio between the minimum and maximum choice of lambda. Default is \code{NULL}, where the ratio is chosen as 1e-2.
#' @param ncov.lambda.weight Weight of the penalty lambda applied to the first \code{ncov} covariates. Default is 0 such that the first \code{ncov} covariates are not penalized.
#' @param a A scalar between 0 and 1: \code{a} is the weight for lasso penalty while \code{1-a} is the weight for ridge penalty.
#' @param mu Value of penalty for the augmented Lagrangian
#' @param maxiter Number of iterations needed for the outer loop of the augmented Lagrangian algorithm.
#' @param ncv Folds of cross-validation. Use \code{NULL} if cross-validation is not wanted.
#' @param ncore Number of cores for parallel computing for cross-validation. Default is 1.
#' @param intercept \code{TRUE} or \code{FALSE}, indicating whether an intercept should be estimated.
#' @param foldid A vector of fold indicator. Default is \code{NULL}.
#' @param step2 \code{TRUE} or \code{FALSE}, indicating whether a second-stage feature selection for specific ratios should be performed for the features selected by the main lasso algorithm. Will only be performed if cross validation is enabled.
#' @param progress \code{TRUE} or \code{FALSE}, indicating whether printing progress bar as the algorithm runs.
#' @param plot \code{TRUE} or \code{FALSE}, indicating whether returning plots of model fitting.
#' @return A list with path-specific estimates (beta), path (lambda), and others. Details can be found in \code{README.md}.
#' @author Teng Fei. Email: feit1@mskcc.org
#' @references Fei T, Funnell T, Waters N, Raj SS et al. Scalable Log-ratio Lasso Regression Enhances Microbiome Feature Selection for Predictive Models. bioRxiv 2023.05.02.538599.
#' 
#' @examples 
#' 
#' set.seed(23420)
#' 
#' # Continuous outcome
#' dat <- simu(n=50,p=30,model="linear")
#' fit <- FLORAL(dat$xcount,dat$y,family="gaussian",ncv=2,progress=FALSE,step2=TRUE)
#' 
#' # Binary outcome
#' # dat <- simu(n=50,p=30,model="binomial")
#' # fit <- FLORAL(dat$xcount,dat$y,family="binomial",progress=FALSE,step2=TRUE)
#' 
#' # Survival outcome
#' # dat <- simu(n=50,p=30,model="cox")
#' # fit <- FLORAL(dat$xcount,survival::Surv(dat$t,dat$d),family="cox",progress=FALSE,step2=TRUE)
#' 
#' # Competing risks outcome
#' # dat <- simu(n=50,p=30,model="finegray")
#' # fit <- FLORAL(dat$xcount,survival::Surv(dat$t,dat$d,type="mstate"),failcode=1,
#' #               family="finegray",progress=FALSE,step2=FALSE)
#' 
#' @import Rcpp ggplot2 survival glmnet dplyr doParallel foreach doRNG parallel
#' @importFrom survcomp concordance.index
#' @importFrom reshape melt
#' @importFrom utils combn
#' @importFrom grDevices rainbow
#' @importFrom caret createFolds
#' @importFrom stats dist rbinom rexp rmultinom rnorm runif sd step glm binomial gaussian na.omit
#' @useDynLib FLORAL
#' @export

FLORAL <- function(x,
                   y,
                   ncov=0,
                   family="gaussian",
                   longitudinal=FALSE,
                   id=NULL,
                   tobs=NULL,
                   failcode=NULL,
                   corstr="exchangeable",
                   scalefix=FALSE,
                   scalevalue=1,
                   pseudo=1,
                   length.lambda=100,
                   lambda.min.ratio=NULL,
                   ncov.lambda.weight=0,
                   a=1,
                   mu=1,
                   maxiter=100,
                   ncv=5,
                   ncore=1,
                   intercept=FALSE,
                   foldid=NULL,
                   step2=TRUE,
                   progress=TRUE,
                   plot=TRUE){
  
  if (ncov < 0){
    stop("`ncov` must be a non-negative integer.")
  }
  
  if (a < 0){
    stop("`a` must be non-negative.")
  }else if (a <= 1){
    if (progress) cat(paste0("Using elastic net with a=",a,"."))
  }else if (a > 1){
    if (longitudinal & family %in% c("gaussian","binomial")){
      if (progress) cat(paste0("Using SCAD with a=",a,"."))
    }else{
      stop("`a`>1 is not yet supported for non-GEE models.")
    }
  }
  
  if (!is.null(pseudo)){
    
    x[,(ncov+1):ncol(x)] <- log(x[,(ncov+1):ncol(x)]+pseudo)
    
  }
  
  if (is.null(colnames(x))){
    colnames(x) = 1:ncol(x)
  }
  
  if (family == "gaussian"){
    
    if (longitudinal){
      
      if (is.null(id)){
        
        stop("`id` must be specified if `longitudinal` is TRUE and `family` is gaussian.")
        
      }else{
        
        if (intercept){
          
          x0 <- cbind(1,x)
          ncov <- ncov + 1
          if (!is.null(colnames(x))){
            colnames(x0) <- c("Intercept",colnames(x))
          }
          x <- x0
        }
        
        res <- LogRatioGEE(x,
                           y,
                           id,
                           ncov,
                           intercept,
                           family,
                           corstr,
                           scalefix,
                           scalevalue,
                           length.lambda,
                           lambda.min.ratio,
                           ncov.lambda.weight,
                           a,
                           mu,
                           ncv,
                           foldid,
                           step2,
                           progress,
                           plot,
                           ncore)
        
      }
      
    }else{
      res <- LogRatioLasso(x,
                           y,
                           ncov,
                           length.lambda,
                           lambda.min.ratio,
                           ncov.lambda.weight,
                           a,
                           mu,
                           maxiter,
                           ncv,
                           intercept,
                           foldid,
                           step2,
                           progress,
                           plot,
                           ncore=ncore)
    }
    
  }else if(family == "binomial"){
    
    if (longitudinal){
      
      if (is.null(id)){
        
        stop("`id` must be specified if `longitudinal` is TRUE and `family` is binomial.")
        
      }else{
        
        if (intercept){
          
          x0 <- cbind(1,x)
          ncov <- ncov + 1
          if (!is.null(colnames(x))){
            colnames(x0) <- c("Intercept",colnames(x))
          }
          x <- x0
        }
        
        res <- LogRatioGEE(x,
                           y,
                           id,
                           ncov,
                           intercept,
                           family,
                           corstr,
                           scalefix,
                           scalevalue,
                           length.lambda,
                           lambda.min.ratio,
                           ncov.lambda.weight,
                           a,
                           mu,
                           ncv,
                           foldid,
                           step2,
                           progress,
                           plot,
                           ncore)
        
      }
      
    }else{
      res <- LogRatioLogisticLasso(x,
                                   y,
                                   ncov,
                                   length.lambda,
                                   lambda.min.ratio,
                                   ncov.lambda.weight,
                                   a,
                                   mu,
                                   maxiter,
                                   ncv,
                                   foldid,
                                   step2,
                                   progress,
                                   plot,
                                   ncore=ncore)
    }
    
  }else if(family == "cox"){
    
    if (longitudinal){
      
      if (is.null(id) | is.null(tobs)){
        
        stop("`id` and `tobs` must be specified if `longitudinal` is TRUE and `family` is cox.")
        
      }else{
        
        if (min(tobs) < 0){
          warning("Negative `tobs` is detected. Please consider adjusting the time origin.")
        }
        
        xidt <- data.frame(cbind(x,id,tobs))
        yid <- data.frame(t=y[,1],d=y[,2],id=unique(id))
        xy <- xidt %>% left_join(yid,by="id") %>% group_by(id) %>% 
          mutate(tstart = tobs,
                 tstop = ifelse(row_number()==n(),t,lead(tobs)),
                 dd = ifelse(row_number()==n(),.data$d,0)) %>% 
          filter(.data$tstart < .data$tstop)
        
        newx <- as.matrix(xy[,colnames(x)])
        newy <- Surv(xy$tstart,xy$tstop,xy$dd)
        newid <- xy$id
        
        res <- LogRatioTDCoxLasso(newx,
                                  newy,
                                  newid,
                                  ncov,
                                  length.lambda,
                                  lambda.min.ratio,
                                  ncov.lambda.weight,
                                  a,
                                  mu,
                                  maxiter,
                                  ncv,
                                  foldid,
                                  step2,
                                  progress,
                                  plot,
                                  ncore=ncore)
        
      }
      
    }else{
      
      res <- LogRatioCoxLasso(x,
                              y,
                              ncov,
                              length.lambda,
                              lambda.min.ratio,
                              ncov.lambda.weight,
                              a,
                              mu,
                              maxiter,
                              ncv,
                              foldid,
                              step2,
                              progress,
                              plot,
                              ncore=ncore)
      
    }
  }else if (family == "finegray"){
    
    if (longitudinal){
      
      if (is.null(id) | is.null(tobs)){
        
        stop("`id` and `tobs` must be specified if `longitudinal` is TRUE and `family` is finegray.")
        
      }else{
        
        if (min(tobs) < 0){
          warning("Negative `tobs` is detected. Please consider adjusting the time origin.")
        }
        
        xidt <- data.frame(cbind(x,id,tobs))
        yid <- data.frame(t=y[,1],d=y[,2],id=unique(id))
        xy <- xidt %>% left_join(yid,by="id") %>% group_by(id) %>% 
          mutate(tstart = tobs,
                 tstop = ifelse(row_number()==n(),t,lead(tobs)),
                 dd = ifelse(row_number()==n(),.data$d,0)) %>% 
          filter(.data$tstart < .data$tstop)
        
        if (is.null(failcode)){
          warning("`failcode` is `NULL`. Using the first failure type as default")
          failcode = 1
          df_FG <- finegray(Surv(tstart, tstop, dd, type="mstate") ~ ., id=id, data=xy, etype=failcode,timefix = FALSE)
        }else{
          df_FG <- finegray(Surv(tstart, tstop, dd, type="mstate") ~ ., id=id, data=xy, etype=failcode,timefix = FALSE)
        }
        
        newx <- as.matrix(df_FG[,colnames(x)])
        newy <- Surv(df_FG$fgstart,df_FG$fgstop,df_FG$fgstatus)
        weight <- df_FG$fgwt
        newid <- df_FG$id
        
        res <- LogRatioFGLasso(newx,
                               newy,
                               newid,
                               weight,
                               ncov,
                               length.lambda,
                               lambda.min.ratio,
                               ncov.lambda.weight,
                               a,
                               mu,
                               maxiter,
                               ncv,
                               foldid,
                               step2,
                               progress,
                               plot,
                               ncore=ncore)
        
      }
      
    }else{
      
      xy <- data.frame(cbind(x,y))
      xy$id <- 1:nrow(xy)
      
      if (is.null(failcode)){
        warning("`failcode` is `NULL`. Using the first failure type as default")
        failcode = 1
        df_FG <- finegray(Surv(time,status,type="mstate") ~ ., data=xy, etype=failcode,timefix = FALSE)
      }else{
        df_FG <- finegray(Surv(time,status,type="mstate") ~ ., data=xy, etype=failcode,timefix = FALSE)
      }
      
      newx <- as.matrix(df_FG[,colnames(x)])
      newy <- Surv(df_FG$fgstart,df_FG$fgstop,df_FG$fgstatus)
      weight <- df_FG$fgwt
      newid <- df_FG$id
      
      res <- LogRatioFGLasso(newx,
                             newy,
                             newid,
                             weight,
                             ncov,
                             length.lambda,
                             lambda.min.ratio,
                             ncov.lambda.weight,
                             a,
                             mu,
                             maxiter,
                             ncv,
                             foldid,
                             step2,
                             progress,
                             plot,
                             ncore=ncore)
      
    }
    
  }
  
  res$loss <- NULL
  res$mse <- NULL
  res$best.beta$min <- res$best.beta$min.mse
  res$best.beta$`1se` <- res$best.beta$add.1se
  res$best.beta$min.mse <- NULL
  res$best.beta$add.1se <- NULL
  
  if (step2){
    
    if (longitudinal & family %in% c("gaussian","binomial")){ # For GEE
      
      res$selected.feature <- list(min=names(res$best.beta$min)[which(res$best.beta$min!=0)],
                                   `1se`=names(res$best.beta$`1se`)[which(res$best.beta$`1se`!=0)],
                                   min.2stage=sort(unique(as.vector(res$step2.feature.min))),
                                   `1se.2stage`=sort(unique(as.vector(res$step2.feature.1se))))
      
      res$step2.ratios <- list(min=character(0),
                               `1se`=character(0),
                               min.name=character(0),
                               `1se.name`=character(0))
      
      res$step2.tables <- list(min=character(0),
                               `1se`=character(0))
      
      
      if (length(res$selected.feature$min.2stage)>0){
        
        namemat <- matrix(res$step2.feature.min,nrow=2)
        res$step2.ratios$min <- as.vector(na.omit(apply(namemat,2,function(x) ifelse(sum(is.na(x))==0,paste(x,collapse ="/"),NA))))
        res$step2.ratios$min.name <- res$step2.feature.min
        res$step2.tables$min <- res$step2fit.min
        
      }
      
      if (length(res$selected.feature$`1se.2stage`)>0){
        
        namemat <- matrix(res$step2.feature.1se,nrow=2)
        res$step2.ratios$`1se` <- as.vector(na.omit(apply(namemat,2,function(x) ifelse(sum(is.na(x))==0,paste(x,collapse ="/"),NA))))
        res$step2.ratios$`1se.name` <- res$step2.feature.1se
        res$step2.tables$`1se` <- res$step2fit.1se
        
      }
      
    }else{ # For other models
      
      
      res$selected.feature <- list(min=names(res$best.beta$min)[which(res$best.beta$min!=0)],
                                   `1se`=names(res$best.beta$`1se`)[which(res$best.beta$`1se`!=0)],
                                   min.2stage=as.vector(na.omit(unique(names(res$best.beta$min)[res$step2.feature.min]))),
                                   `1se.2stage`=as.vector(na.omit(unique(names(res$best.beta$min)[res$step2.feature.1se])))
      )
      
      res$step2.ratios <- list(min=character(0),
                               `1se`=character(0),
                               min.idx=character(0),
                               `1se.idx`=character(0))
      
      res$step2.tables <- list(min=character(0),
                               `1se`=character(0))
      
      if (length(res$selected.feature$min.2stage)>0){
        namemat <- matrix(names(res$best.beta$min)[res$step2.feature.min],nrow=2)
        res$step2.ratios$min <- as.vector(na.omit(apply(namemat,2,function(x) ifelse(sum(is.na(x))==0,paste(x,collapse ="/"),NA))))
        res$step2.ratios$min.idx <- res$step2.feature.min#[,!is.na(colSums(res$step2.feature.min))]
        
        res$step2.tables$`min` <- summary(res$step2fit.min)$coefficients
        
        if ("(Intercept)" %in% rownames(res$step2.tables$min)){
          rownames(res$step2.tables$`min`)[2:(2+length(res$step2.ratios$min)-1)] <- res$step2.ratios$`min`
        }else{
          rownames(res$step2.tables$`min`)[1:length(res$step2.ratios$min)] <- res$step2.ratios$`min`
        }
        
        # res$step2.tables$`min` <- res$step2.tables$`min`[rownames(res$step2.tables$`min`) != "(Intercept)", ]
        
      }
      if (length(res$selected.feature$`1se.2stage`)>0){
        namemat <- matrix(names(res$best.beta$`1se`)[res$step2.feature.1se],nrow=2)
        res$step2.ratios$`1se` <- as.vector(na.omit(apply(namemat,2,function(x) ifelse(sum(is.na(x))==0,paste(x,collapse ="/"),NA))))
        res$step2.ratios$`1se.idx` <- res$step2.feature.1se#[,!is.na(colSums(res$step2.feature.1se))]
        
        res$step2.tables$`1se` <- summary(res$step2fit.1se)$coefficients
        # res$step2.tables$`1se` <- res$step2.tables$`1se`[rownames(res$step2.tables$`1se`) != "(Intercept)", ]
        # rownames(res$step2.tables$`1se`)[colSums(is.na(namemat)) == 0] <- res$step2.ratios$`1se`
        if ("(Intercept)" %in% rownames(res$step2.tables$`1se`)){
          rownames(res$step2.tables$`1se`)[2:(2+length(res$step2.ratios$`1se`)-1)] <- res$step2.ratios$`1se`
        }else{
          rownames(res$step2.tables$`1se`)[1:length(res$step2.ratios$`1se`)] <- res$step2.ratios$`1se`
        }
      }
    }
    
  }else{
    res$selected.feature <- list(min=names(res$best.beta$min)[which(res$best.beta$min!=0)],
                                 `1se`=names(res$best.beta$`1se`)[which(res$best.beta$`1se`!=0)]
    )
  }
  
  res$step2.feature.min <- NULL
  res$step2.feature.1se <- NULL
  
  res$selected.feature <- lapply(res$selected.feature,sort)
  
  return(res)
  
}

#' Summarizing selected compositional features over multiple cross validations
#' 
#' @description Summarizing \code{FLORAL} outputs from multiple random k-fold cross validations
#' @param mcv Number of random `ncv`-fold cross-validation to be performed.
#' @param ncore Number of cores used for parallel computation. Default is to use only 1 core.
#' @param seed A random seed for reproducibility of the results. By default the seed is the numeric form of \code{Sys.Date()}.
#' @param x Feature matrix, where rows specify subjects and columns specify features. The first \code{ncov} columns should be patient characteristics and the rest columns are microbiome absolute counts corresponding to various taxa. If \code{x} contains longitudinal data, the rows must be sorted in the same order of the subject IDs used in \code{y}.
#' @param y Outcome. For a continuous or binary outcome, \code{y} is a vector. For survival outcome, \code{y} is a \code{Surv} object.
#' @param ncov An integer indicating the number of first \code{ncov} columns in \code{x} that will not be subject to the zero-sum constraint.
#' @param family Available options are \code{gaussian}, \code{binomial}, \code{cox}, \code{finegray}.
#' @param longitudinal \code{TRUE} or \code{FALSE}, indicating whether longitudinal data matrix is specified for input \code{x}. (\code{Longitudinal=TRUE} and \code{family="cox"} or \code{"finegray"} will fit a time-dependent covariate model. \code{Longitudinal=TRUE} and \code{family="gaussian"} or \code{"binomial"} will fit a GEE model.)
#' @param id If \code{longitudinal} is \code{TRUE}, \code{id} specifies subject IDs corresponding to the rows of input \code{x}.
#' @param tobs If \code{longitudinal} is \code{TRUE}, \code{tobs} specifies time points corresponding to the rows of input \code{x}.
#' @param failcode If \code{family = finegray}, \code{failcode} specifies the failure type of interest. This must be a positive integer.
#' @param corstr If a GEE model is specified, then \code{corstr} is the corresponding working correlation structure. Options are \code{independence}, \code{exchangeable}, \code{AR-1} and \code{unstructured}.
#' @param scalefix \code{TRUE} or \code{FALSE}, indicating whether the scale parameter is estimated or fixed if a GEE model is specified.
#' @param scalevalue Specify the scale parameter if \code{scalefix=TRUE}.
#' @param pseudo Pseudo count to be added to \code{x} before taking log-transformation
#' @param length.lambda Number of penalty parameters used in the path
#' @param lambda.min.ratio Ratio between the minimum and maximum choice of lambda. Default is \code{NULL}, where the ratio is chosen as 1e-2.
#' @param ncov.lambda.weight Weight of the penalty lambda applied to the first \code{ncov} covariates. Default is 0 such that the first \code{ncov} covariates are not penalized.
#' @param a A scalar between 0 and 1: \code{a} is the weight for lasso penalty while \code{1-a} is the weight for ridge penalty.
#' @param mu Value of penalty for the augmented Lagrangian
#' @param maxiter Number of iterations needed for the outer loop of the augmented Lagrangian algorithm.
#' @param ncv Folds of cross-validation. Use \code{NULL} if cross-validation is not wanted.
#' @param intercept \code{TRUE} or \code{FALSE}, indicating whether an intercept should be estimated.
#' @param step2 \code{TRUE} or \code{FALSE}, indicating whether a second-stage feature selection for specific ratios should be performed for the features selected by the main lasso algorithm. Will only be performed if cross validation is enabled.
#' @param progress \code{TRUE} or \code{FALSE}, indicating whether printing progress bar as the algorithm runs.
#' @param plot \code{TRUE} or \code{FALSE}, indicating whether returning summary plots of selection probability for taxa features.
#' @return A list with relative frequencies of a certain feature being selected over \code{mcv} \code{ncv}-fold cross-validations.
#' @author Teng Fei. Email: feit1@mskcc.org
#' @references Fei T, Funnell T, Waters N, Raj SS et al. Scalable Log-ratio Lasso Regression Enhances Microbiome Feature Selection for Predictive Models. bioRxiv 2023.05.02.538599.
#' 
#' @examples 
#' 
#' set.seed(23420)
#' 
#' dat <- simu(n=50,p=30,model="linear")
#' fit <- mcv.FLORAL(mcv=2,ncore=1,x=dat$xcount,y=dat$y,ncv=2,progress=FALSE,step2=TRUE,plot=FALSE)
#' 
#' @import Rcpp ggplot2 survival glmnet dplyr doParallel foreach parallel
#' @importFrom survcomp concordance.index
#' @importFrom reshape melt
#' @importFrom utils combn
#' @importFrom grDevices rainbow
#' @importFrom caret createFolds
#' @importFrom stats dist rbinom rexp rmultinom rnorm runif sd step glm binomial gaussian na.omit
#' @useDynLib FLORAL
#' @export

mcv.FLORAL <- function(mcv=10,
                       ncore=1,
                       seed=NULL,
                       x,
                       y,
                       ncov=0,
                       family="gaussian",
                       longitudinal=FALSE,
                       id=NULL,
                       tobs=NULL,
                       failcode=NULL,
                       corstr="exchangeable",
                       scalefix=FALSE,
                       scalevalue=1,
                       pseudo=1,
                       length.lambda=100,
                       lambda.min.ratio=NULL,
                       ncov.lambda.weight=0,
                       a=1,
                       mu=1,
                       maxiter=100,
                       ncv=5,
                       intercept=FALSE,
                       step2=TRUE,
                       progress=TRUE,
                       plot=TRUE){
  
  if (is.null(ncv) | ncv < 2) stop("Number of folds `ncv` must be larger than one for cross validation.")
  
  if (mcv <= 1){
    
    stop("`mcv` is less than 2. Please directly use `FLORAL` function.")
    
  }else if (mcv >= 2){
    
    if (ncore == 1){
      
      if (progress) warning("Using 1 core for computation.")
      
      FLORAL.res <- list()
      if (is.null(seed)) seed <- as.numeric(Sys.Date())
      set.seed(seed)
      
      for (i in 1:mcv){
        
        if (progress) cat(paste0("Random ",ncv,"-fold cross-validation: ",i,"\n"))
        
        FLORAL.res[[i]] <- FLORAL(x,
                                  y,
                                  ncov,
                                  family,
                                  longitudinal,
                                  id,
                                  tobs,
                                  failcode,
                                  corstr,
                                  scalefix,
                                  scalevalue,
                                  pseudo,
                                  length.lambda,
                                  lambda.min.ratio,
                                  ncov.lambda.weight,
                                  a,
                                  mu,
                                  maxiter,
                                  ncv,
                                  ncore=1,
                                  intercept,
                                  foldid=NULL,
                                  step2,
                                  progress=FALSE,
                                  plot=FALSE)
        
      }
      
      if (longitudinal & family %in% c("gaussian","binomial")){
        
        res <- list(min=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$min)))/mcv,
                    `1se`=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$`1se`)))/mcv,
                    min.2stage=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$min.2stage)))/mcv,
                    `1se.2stage`=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$`1se.2stage`)))/mcv,
                    min.2stage.ratios=table(unlist(lapply(FLORAL.res,function(x) names(x$step2.tables$min))))/mcv,
                    `1se.2stage.ratios`=table(unlist(lapply(FLORAL.res,function(x) names(x$step2.tables$`1se`))))/mcv,
                    mcv=mcv,
                    seed=seed)
        
        res$min.coef = Reduce(`+`,lapply(FLORAL.res,function(x) x$best.beta$min))[names(res$min)]/Reduce(`+`,lapply(lapply(FLORAL.res, function(x) x$best.beta$min), function(x) x != 0))[names(res$min)]
        res$`1se.coef` = Reduce(`+`,lapply(FLORAL.res,function(x) x$best.beta$`1se`))[names(res$`1se`)]/Reduce(`+`,lapply(lapply(FLORAL.res, function(x) x$best.beta$`1se`), function(x) x != 0))[names(res$`1se`)]
        res$min.2stage.ratios.coef = colMeans(bind_rows(lapply(FLORAL.res,function(x) if(length(x$step2.tables$min)>0) x$step2.tables$min)),na.rm=TRUE)[names(res$min.2stage.ratios)]
        res$`1se.2stage.ratios.coef` = colMeans(bind_rows(lapply(FLORAL.res,function(x) if(length(x$step2.tables$`1se`)>0) x$step2.tables$`1se`)),na.rm=TRUE)[names(res$`1se.2stage.ratios`)]
        
        
      }else{
        
        res <- list(min=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$min)))/mcv,
                    `1se`=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$`1se`)))/mcv,
                    min.2stage=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$min.2stage)))/mcv,
                    `1se.2stage`=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$`1se.2stage`)))/mcv,
                    min.2stage.ratios=table(unlist(lapply(FLORAL.res,function(x) rownames(x$step2.tables$min))))/mcv,
                    `1se.2stage.ratios`=table(unlist(lapply(FLORAL.res,function(x) rownames(x$step2.tables$`1se`))))/mcv,
                    mcv=mcv,
                    seed=seed)
        
        res$min.coef = Reduce(`+`,lapply(FLORAL.res,function(x) x$best.beta$min))[names(res$min)]/Reduce(`+`,lapply(lapply(FLORAL.res, function(x) x$best.beta$min), function(x) x != 0))[names(res$min)]
        res$`1se.coef` = Reduce(`+`,lapply(FLORAL.res,function(x) x$best.beta$`1se`))[names(res$`1se`)]/Reduce(`+`,lapply(lapply(FLORAL.res, function(x) x$best.beta$`1se`), function(x) x != 0))[names(res$`1se`)]
        res$min.2stage.ratios.coef = colMeans(bind_rows(lapply(FLORAL.res,function(x) if(length(x$step2.tables$min)>0) x$step2.tables$min[,1])),na.rm=TRUE)[names(res$min.2stage.ratios)]
        res$`1se.2stage.ratios.coef` = colMeans(bind_rows(lapply(FLORAL.res,function(x) if(length(x$step2.tables$`1se`)>0) x$step2.tables$`1se`[,1])),na.rm=TRUE)[names(res$`1se.2stage.ratios`)]
        
        if (family %in% c("cox","finegray")){
          res$min.2stage.ratios.p = colMeans(bind_rows(lapply(FLORAL.res,function(x) if(length(x$step2.tables$min)>0) x$step2.tables$min[,5])),na.rm=TRUE)[names(res$min.2stage.ratios)]
          res$`1se.2stage.ratios.p` = colMeans(bind_rows(lapply(FLORAL.res,function(x) if(length(x$step2.tables$`1se`)>0) x$step2.tables$`1se`[,5])),na.rm=TRUE)[names(res$`1se.2stage.ratios`)]
        }else if (family %in% c("gaussian","binomial")){
          res$min.2stage.ratios.p = colMeans(bind_rows(lapply(FLORAL.res,function(x) if(length(x$step2.tables$min)>0) x$step2.tables$min[,4])),na.rm=TRUE)[names(res$min.2stage.ratios)]
          res$`1se.2stage.ratios.p` = colMeans(bind_rows(lapply(FLORAL.res,function(x) if(length(x$step2.tables$`1se`)>0) x$step2.tables$`1se`[,4])),na.rm=TRUE)[names(res$`1se.2stage.ratios`)]
        }
        
      }
      
    }else if (ncore > 1){
      
      if (progress) warning(paste0("Using ", ncore ," core for computation."))
      
      cl <- makeCluster(ncore)
      registerDoParallel(cl)
      if (is.null(seed)) seed <- as.numeric(Sys.Date())
      registerDoRNG(seed=seed)
      
      FLORAL.res <- foreach(i=1:mcv) %dopar% {
        
        FLORAL(x,
               y,
               ncov,
               family,
               longitudinal,
               id,
               tobs,
               failcode,
               corstr,
               scalefix,
               scalevalue,
               pseudo,
               length.lambda,
               lambda.min.ratio,
               ncov.lambda.weight,
               a,
               mu,
               maxiter,
               ncv,
               ncore=1,
               intercept,
               foldid=NULL,
               step2,
               progress=FALSE,
               plot=FALSE)
        
      }
      
      stopCluster(cl)
      
      if (longitudinal & family %in% c("gaussian","binomial")){
        
        res <- list(min=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$min)))/mcv,
                    `1se`=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$`1se`)))/mcv,
                    min.2stage=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$min.2stage)))/mcv,
                    `1se.2stage`=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$`1se.2stage`)))/mcv,
                    min.2stage.ratios=table(unlist(lapply(FLORAL.res,function(x) names(x$step2.tables$min))))/mcv,
                    `1se.2stage.ratios`=table(unlist(lapply(FLORAL.res,function(x) names(x$step2.tables$`1se`))))/mcv,
                    mcv=mcv,
                    seed=seed)
        
        res$min.coef = Reduce(`+`,lapply(FLORAL.res,function(x) x$best.beta$min))[names(res$min)]/Reduce(`+`,lapply(lapply(FLORAL.res, function(x) x$best.beta$min), function(x) x != 0))[names(res$min)]
        res$`1se.coef` = Reduce(`+`,lapply(FLORAL.res,function(x) x$best.beta$`1se`))[names(res$`1se`)]/Reduce(`+`,lapply(lapply(FLORAL.res, function(x) x$best.beta$`1se`), function(x) x != 0))[names(res$`1se`)]
        res$min.2stage.ratios.coef = colMeans(bind_rows(lapply(FLORAL.res,function(x) if(length(x$step2.tables$min)>0) x$step2.tables$min)),na.rm=TRUE)[names(res$min.2stage.ratios)]
        res$`1se.2stage.ratios.coef` = colMeans(bind_rows(lapply(FLORAL.res,function(x) if(length(x$step2.tables$`1se`)>0) x$step2.tables$`1se`)),na.rm=TRUE)[names(res$`1se.2stage.ratios`)]
        
        
      }else{
        
        res <- list(min=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$min)))/mcv,
                    `1se`=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$`1se`)))/mcv,
                    min.2stage=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$min.2stage)))/mcv,
                    `1se.2stage`=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$`1se.2stage`)))/mcv,
                    min.2stage.ratios=table(unlist(lapply(FLORAL.res,function(x) rownames(x$step2.tables$min))))/mcv,
                    `1se.2stage.ratios`=table(unlist(lapply(FLORAL.res,function(x) rownames(x$step2.tables$`1se`))))/mcv,
                    mcv=mcv,
                    seed=seed)
        
        res$min.coef = Reduce(`+`,lapply(FLORAL.res,function(x) x$best.beta$min))[names(res$min)]/Reduce(`+`,lapply(lapply(FLORAL.res, function(x) x$best.beta$min), function(x) x != 0))[names(res$min)]
        res$`1se.coef` = Reduce(`+`,lapply(FLORAL.res,function(x) x$best.beta$`1se`))[names(res$`1se`)]/Reduce(`+`,lapply(lapply(FLORAL.res, function(x) x$best.beta$`1se`), function(x) x != 0))[names(res$`1se`)]
        res$min.2stage.ratios.coef = colMeans(bind_rows(lapply(FLORAL.res,function(x) if(length(x$step2.tables$min)>0) x$step2.tables$min[,1])),na.rm=TRUE)[names(res$min.2stage.ratios)]
        res$`1se.2stage.ratios.coef` = colMeans(bind_rows(lapply(FLORAL.res,function(x) if(length(x$step2.tables$`1se`)>0) x$step2.tables$`1se`[,1])),na.rm=TRUE)[names(res$`1se.2stage.ratios`)]
        
        if (family %in% c("cox","finegray")){
          res$min.2stage.ratios.p = colMeans(bind_rows(lapply(FLORAL.res,function(x) if(length(x$step2.tables$min)>0) x$step2.tables$min[,5])),na.rm=TRUE)[names(res$min.2stage.ratios)]
          res$`1se.2stage.ratios.p` = colMeans(bind_rows(lapply(FLORAL.res,function(x) if(length(x$step2.tables$`1se`)>0) x$step2.tables$`1se`[,5])),na.rm=TRUE)[names(res$`1se.2stage.ratios`)]
        }else if (family %in% c("gaussian","binomial")){
          res$min.2stage.ratios.p = colMeans(bind_rows(lapply(FLORAL.res,function(x) if(length(x$step2.tables$min)>0) x$step2.tables$min[,4])),na.rm=TRUE)[names(res$min.2stage.ratios)]
          res$`1se.2stage.ratios.p` = colMeans(bind_rows(lapply(FLORAL.res,function(x) if(length(x$step2.tables$`1se`)>0) x$step2.tables$`1se`[,4])),na.rm=TRUE)[names(res$`1se.2stage.ratios`)]
        }
        
      }
      
    }
    
  }
  
  if (plot){
    
    df_plot <- data.frame(taxa=names(sort(res$min.2stage)),
                          prob=as.vector(sort(res$min.2stage)))
    df_plot$Avg.coef <- res$min.coef[df_plot$taxa]
    df_plot$coefsign <- sign(df_plot$Avg.coef)
    
    res$p_min <- ggplot(df_plot, aes(y=.data$taxa,fill=.data$Avg.coef)) + 
      geom_bar(aes(weight=.data$prob),color="darkgrey") +
      scale_y_discrete(limits = df_plot$taxa[order(df_plot$prob,decreasing = T)]) +
      xlab("Probability of being selected") +
      ylab("Taxa") +
      ggtitle(expression(paste(lambda," = ",lambda["min"]))) +
      xlim(0,1) +
      scale_fill_gradient2(low="blue",high="red")+
      theme_bw()
    
    df_plot <- data.frame(taxa=names(sort(res$`1se.2stage`)),
                          prob=as.vector(sort(res$`1se.2stage`)))
    df_plot$Avg.coef <- res$`1se.coef`[df_plot$taxa]
    df_plot$coefsign <- sign(df_plot$Avg.coef)
    
    res$p_1se <- ggplot(df_plot, aes(y=.data$taxa,fill=.data$Avg.coef)) + 
      geom_bar(aes(weight=.data$prob),color="darkgrey") +
      scale_y_discrete(limits = df_plot$taxa[order(df_plot$prob,decreasing = T)]) +
      xlab("Probability of being selected") +
      ylab("Taxa") +
      ggtitle(expression(paste(lambda," = ",lambda["1se"]))) +
      xlim(0,1) +
      scale_fill_gradient2(low="blue",high="red")+
      theme_bw()
    
    df_plot <- data.frame(taxa=names(sort(res$min.2stage.ratios)),
                          prob=as.vector(sort(res$min.2stage.ratios)))
    df_plot$Avg.coef <- res$min.2stage.ratios.coef[df_plot$taxa]
    df_plot$coefsign <- sign(df_plot$Avg.coef)
    
    if (longitudinal & family %in% c("gaussian","binomial")){
      
      res$p_min_ratio <- ggplot(df_plot, aes(y=.data$taxa,fill=.data$Avg.coef)) + 
        geom_bar(aes(weight=.data$prob),color="darkgrey") +
        # geom_point(aes(x=.data$prob,y=.data$taxa,size=.data$p),color="black",alpha=0.5) +
        scale_shape_binned() +
        scale_y_discrete(limits = df_plot$taxa[order(df_plot$prob,decreasing = T)]) +
        xlab("Probability of being selected") +
        ylab("2 stage model covariates") +
        ggtitle(expression(paste(lambda," = ",lambda["min"]))) +
        xlim(0,1) +
        scale_fill_gradient2(low="blue",high="red")+
        theme_bw()
      
    }else{
      
      df_plot$p <- -log10(res$min.2stage.ratios.p[df_plot$taxa])
      
      res$p_min_ratio <- ggplot(df_plot, aes(y=.data$taxa,fill=.data$Avg.coef)) + 
        geom_bar(aes(weight=.data$prob),color="darkgrey") +
        geom_point(aes(x=.data$prob,y=.data$taxa,size=.data$p),color="black",alpha=0.5) +
        scale_shape_binned() +
        scale_y_discrete(limits = df_plot$taxa[order(df_plot$prob,decreasing = T)]) +
        xlab("Probability of being selected") +
        ylab("2 stage model covariates") +
        ggtitle(expression(paste(lambda," = ",lambda["min"]))) +
        xlim(0,1) +
        scale_fill_gradient2(low="blue",high="red")+
        theme_bw() + 
        guides(size=guide_legend(title="Avg.-log10(p)"))
      
    }
    
    
    df_plot <- data.frame(taxa=names(sort(res$`1se.2stage.ratios`)),
                          prob=as.vector(sort(res$`1se.2stage.ratios`)))
    df_plot$Avg.coef <- res$`1se.2stage.ratios.coef`[df_plot$taxa]
    df_plot$coefsign <- sign(df_plot$Avg.coef)
    
    if (longitudinal & family %in% c("gaussian","binomial")){
      
      res$p_1se_ratio <- ggplot(df_plot, aes(y=.data$taxa,fill=.data$Avg.coef)) + 
        geom_bar(aes(weight=.data$prob),color="darkgrey") +
        # geom_point(aes(x=.data$prob,y=.data$taxa,size=.data$p),color="black",alpha=0.5) +
        scale_y_discrete(limits = df_plot$taxa[order(df_plot$prob,decreasing = T)]) +
        xlab("Probability of being selected") +
        ylab("2 stage model covariates") +
        ggtitle(expression(paste(lambda," = ",lambda["1se"]))) +
        xlim(0,1) +
        scale_fill_gradient2(low="blue",high="red")+
        theme_bw() 
      
    }else{
      
      df_plot$p <- -log10(res$`1se.2stage.ratios.p`[df_plot$taxa])
      
      res$p_1se_ratio <- ggplot(df_plot, aes(y=.data$taxa,fill=.data$Avg.coef)) + 
        geom_bar(aes(weight=.data$prob),color="darkgrey") +
        geom_point(aes(x=.data$prob,y=.data$taxa,size=.data$p),color="black",alpha=0.5) +
        scale_y_discrete(limits = df_plot$taxa[order(df_plot$prob,decreasing = T)]) +
        xlab("Probability of being selected") +
        ylab("2 stage model covariates") +
        ggtitle(expression(paste(lambda," = ",lambda["1se"]))) +
        xlim(0,1) +
        scale_fill_gradient2(low="blue",high="red")+
        theme_bw() + 
        guides(size=guide_legend(title="Avg.-log10(p)"))
      
    }
    
    
  }
  
  return(res)
  
}


#' Comparing prediction performances under different choices of weights for lasso/ridge penalty
#' 
#' @description Summarizing \code{FLORAL} outputs from various choices of \code{a}
#' @param a vector of scalars between 0 and 1 for comparison.
#' @param ncore Number of cores used for parallel computation. Default is to use only 1 core.
#' @param seed A random seed for reproducibility of the results. By default the seed is the numeric form of \code{Sys.Date()}.
#' @param x Feature matrix, where rows specify subjects and columns specify features. The first \code{ncov} columns should be patient characteristics and the rest columns are microbiome absolute counts corresponding to various taxa. If \code{x} contains longitudinal data, the rows must be sorted in the same order of the subject IDs used in \code{y}.
#' @param y Outcome. For a continuous or binary outcome, \code{y} is a vector. For survival outcome, \code{y} is a \code{Surv} object.
#' @param ncov An integer indicating the number of first \code{ncov} columns in \code{x} that will not be subject to the zero-sum constraint.
#' @param family Available options are \code{gaussian}, \code{binomial}, \code{cox}, \code{finegray}.
#' @param longitudinal \code{TRUE} or \code{FALSE}, indicating whether longitudinal data matrix is specified for input \code{x}. (\code{Longitudinal=TRUE} and \code{family="cox"} or \code{"finegray"} will fit a time-dependent covariate model. \code{Longitudinal=TRUE} and \code{family="gaussian"} or \code{"binomial"} will fit a GEE model.)
#' @param id If \code{longitudinal} is \code{TRUE}, \code{id} specifies subject IDs corresponding to the rows of input \code{x}.
#' @param tobs If \code{longitudinal} is \code{TRUE}, \code{tobs} specifies time points corresponding to the rows of input \code{x}.
#' @param failcode If \code{family = finegray}, \code{failcode} specifies the failure type of interest. This must be a positive integer.
#' @param corstr If a GEE model is specified, then \code{corstr} is the corresponding working correlation structure. Options are \code{independence}, \code{exchangeable}, \code{AR-1} and \code{unstructured}.
#' @param scalefix \code{TRUE} or \code{FALSE}, indicating whether the scale parameter is estimated or fixed if a GEE model is specified.
#' @param scalevalue Specify the scale parameter if \code{scalefix=TRUE}.#' @param pseudo Pseudo count to be added to \code{x} before taking log-transformation
#' @param length.lambda Number of penalty parameters used in the path
#' @param lambda.min.ratio Ratio between the minimum and maximum choice of lambda. Default is \code{NULL}, where the ratio is chosen as 1e-2.
#' @param ncov.lambda.weight Weight of the penalty lambda applied to the first \code{ncov} covariates. Default is 0 such that the first \code{ncov} covariates are not penalized.
#' @param mu Value of penalty for the augmented Lagrangian
#' @param maxiter Number of iterations needed for the outer loop of the augmented Lagrangian algorithm.
#' @param ncv Folds of cross-validation. Use \code{NULL} if cross-validation is not wanted.
#' @param intercept \code{TRUE} or \code{FALSE}, indicating whether an intercept should be estimated.
#' @param step2 \code{TRUE} or \code{FALSE}, indicating whether a second-stage feature selection for specific ratios should be performed for the features selected by the main lasso algorithm. Will only be performed if cross validation is enabled.
#' @param progress \code{TRUE} or \code{FALSE}, indicating whether printing progress bar as the algorithm runs.
#' @return A \code{ggplot2} object of cross-validated prediction metric versus \code{lambda}, stratified by \code{a}. Detailed data can be retrieved from the \code{ggplot2} object itself.
#' @author Teng Fei. Email: feit1@mskcc.org
#' @references Fei T, Funnell T, Waters N, Raj SS et al. Scalable Log-ratio Lasso Regression Enhances Microbiome Feature Selection for Predictive Models. bioRxiv 2023.05.02.538599.
#' 
#' @examples 
#' 
#' set.seed(23420)
#' 
#' dat <- simu(n=50,p=30,model="linear")
#' pmetric <- a.FLORAL(a=c(0.1,1),ncore=1,x=dat$xcount,y=dat$y,family="gaussian",ncv=2,progress=FALSE)
#' 
#' @import Rcpp ggplot2 survival glmnet dplyr doParallel foreach doRNG parallel
#' @importFrom survcomp concordance.index
#' @importFrom reshape melt
#' @importFrom utils combn
#' @importFrom grDevices rainbow
#' @importFrom caret createFolds
#' @importFrom stats dist rbinom rexp rmultinom rnorm runif sd step glm binomial gaussian na.omit
#' @useDynLib FLORAL
#' @export

a.FLORAL <- function(a=c(0.1,0.5,1),
                     ncore=1,
                     seed=NULL,
                     x,
                     y,
                     ncov=0,
                     family="gaussian",
                     longitudinal=FALSE,
                     id=NULL,
                     tobs=NULL,
                     failcode=NULL,
                     corstr="exchangeable",
                     scalefix=FALSE,
                     scalevalue=1,
                     pseudo=1,
                     length.lambda=100,
                     lambda.min.ratio=NULL,
                     ncov.lambda.weight=0,
                     mu=1,
                     maxiter=100,
                     ncv=5,
                     intercept=FALSE,
                     step2=FALSE,
                     progress=TRUE){
  
  if (is.null(ncv) | ncv < 2) stop("Number of folds `ncv` must be larger than one for cross validation.")
  
  if (length(a) <= 1){
    
    stop("Length of `a` is less than 2. Please directly use `FLORAL` function.")
    
  }else if (length(a) >= 2){
    
    if (ncore == 1){
      
      if (progress) warning("Using 1 core for computation.")
      
      FLORAL.res <- list()
      if (is.null(seed)) seed <- as.numeric(Sys.Date())
      set.seed(seed)
      
      for (i in 1:length(a)){
        
        if (progress) cat(paste0("Running for a = ",a[i],"\n"))
        
        if (i == 1){
          foldid=NULL
        }else{
          foldid=fit$foldid
        }
        
        fit <- FLORAL(x,
                      y,
                      ncov,
                      family,
                      longitudinal,
                      id,
                      tobs,
                      failcode,
                      corstr,
                      scalefix,
                      scalevalue,
                      pseudo,
                      length.lambda,
                      lambda.min.ratio,
                      ncov.lambda.weight,
                      a = a[i],
                      mu,
                      maxiter,
                      ncv,
                      ncore=1,
                      intercept,
                      foldid=foldid,
                      step2,
                      progress=FALSE,
                      plot=FALSE)
        
        FLORAL.res[[i]] <- data.frame(a=fit$a,
                                      lambda=as.vector(fit$lambda),
                                      cv.metric=fit[[4]],
                                      cv.metricse=fit[[5]]
        )
        
      }
      
      FLORAL.res <- do.call(rbind,FLORAL.res) %>% 
        group_by(a) %>% 
        mutate(min.metric = min(.data$cv.metric),
               min.lambda = .data$lambda[.data$min.metric==.data$cv.metric]) %>% 
        ungroup()
      
      pmetric <- ggplot(aes(x=log(.data$lambda),y=.data$cv.metric,color=as.factor(.data$a)),data=FLORAL.res) + 
        geom_point(size=0.5) +
        geom_hline(aes(yintercept = .data$min.metric,color=as.factor(.data$a)),
                   alpha = 0.5) +
        geom_vline(aes(xintercept = log(.data$min.lambda),color=as.factor(.data$a)),
                   alpha = 0.5) +
        theme_bw() +
        xlab(expression(paste("log(",lambda,")"))) +
        ylab("Cross-validated metric") + 
        guides(color=guide_legend(title="a"))
      
    }else if (ncore > 1){
      
      if (progress) warning(paste0("Using ", ncore ," core for computation."))
      
      cl <- makeCluster(ncore)
      registerDoParallel(cl)
      if (is.null(seed)) seed <- as.numeric(Sys.Date())
      registerDoRNG(seed=seed)
      
      fit0 <- FLORAL(x,
                     y,
                     ncov,
                     family,
                     longitudinal,
                     id,
                     tobs,
                     failcode,
                     corstr,
                     scalefix,
                     scalevalue,
                     pseudo,
                     length.lambda,
                     lambda.min.ratio,
                     ncov.lambda.weight,
                     a = a[1],
                     mu,
                     maxiter,
                     ncv,
                     ncore=1,
                     intercept,
                     foldid=NULL,
                     step2,
                     progress=FALSE,
                     plot=FALSE)
      
      FLORAL.res0 <- data.frame(a=fit0$a,
                                lambda=as.vector(fit0$lambda),
                                cv.metric=fit0[[4]],
                                cv.metricse=fit0[[5]])
      
      FLORAL.res <- foreach(i=2:length(a)) %dopar% {
        
        fit <- FLORAL(x,
                      y,
                      ncov,
                      family,
                      longitudinal,
                      id,
                      tobs,
                      failcode,
                      corstr,
                      scalefix,
                      scalevalue,
                      pseudo,
                      length.lambda,
                      lambda.min.ratio,
                      ncov.lambda.weight,
                      a = a[i],
                      mu,
                      maxiter,
                      ncv,
                      ncore=1,
                      intercept,
                      foldid=fit0$foldid,
                      step2,
                      progress=FALSE,
                      plot=FALSE)
        
        data.frame(a=fit$a,
                   lambda=as.vector(fit$lambda),
                   cv.metric=fit[[4]],
                   cv.metricse=fit[[5]]
        )
        
      }
      
      stopCluster(cl)
      
      FLORAL.res[[length(a)]] <- FLORAL.res0
      
      FLORAL.res <- do.call(rbind,FLORAL.res) %>% 
        group_by(a) %>% 
        mutate(min.metric = min(.data$cv.metric),
               min.lambda = .data$lambda[.data$min.metric==.data$cv.metric]) %>% 
        ungroup()
      
      pmetric <- ggplot(aes(x=log(.data$lambda),y=.data$cv.metric,color=as.factor(.data$a)),data=FLORAL.res) + 
        geom_point(size=0.5) +
        geom_hline(aes(yintercept = .data$min.metric,color=as.factor(.data$a)),
                   alpha = 0.5) +
        geom_vline(aes(xintercept = log(.data$min.lambda),color=as.factor(.data$a)),
                   alpha = 0.5) +
        theme_bw() +
        xlab(expression(paste("log(",lambda,")"))) +
        ylab("Cross-validated metric") + 
        guides(color=guide_legend(title="a"))
      
    }
    
  }
  
  return(pmetric)
  
}
