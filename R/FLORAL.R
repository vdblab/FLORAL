#' FLORAL: Fit Log-ratio lasso regression for clinical outcomes 
#'
#' @description Conduct log-ratio lasso regression for continuous, binary and survival outcomes. 
#' @param x Count data matrix, where rows specify subjects and columns specify features. If `x` contains longitudinal data, the rows must be sorted in the same order of the subject IDs used in `y`.
#' @param y Outcome. For a continuous or binary outcome, `y` is a vector. For survival outcome, `y` is a `Surv` object.
#' @param family Available options are `gaussian`, `binomial`, `cox`, `finegray`.
#' @param longitudinal TRUE or FALSE, indicating whether longitudinal data matrix is specified for input `x`.
#' @param id If `longitudinal` is TRUE, `id` specifies subject IDs corresponding to the rows of input `x`.
#' @param tobs If `longitudinal` is TRUE, `tobs` specifies time points corresponding to the rows of input `x`.
#' @param failcode If `family = finegray`, `failcode` specifies the failure type of interest. This must be a positive integer.
#' @param length.lambda Number of penalty parameters used in the path
#' @param lambda.min.ratio Ratio between the minimum and maximum choice of lambda. Default is `NULL`, where the ratio is chosen as 1e-2 if n < p and 1e-4 otherwise.
#' @param mu Value of penalty for the augmented Lagrangian
#' @param ncv Number of cross-validation runs. Use `NULL` if cross-validation is not wanted.
#' @param intercept TRUE or FALSE, indicating whether an intercept should be estimated.
#' @param foldid A vector of fold indicator. Default is `NULL`.
#' @param step2 TRUE or FALSE, indicating whether a stepwise feature selection should be performed for features selected by the main lasso algorithm. Will only be performed if cross validation is enabled.
#' @param progress TRUE or FALSE, indicating whether printing progress bar as the algorithm runs.
#' @param plot TRUE or FALSE, indicating whether returning plots of model fitting.
#' @return A list with path-specific estimates (beta), path (lambda), and many others.
#' @author Teng Fei. Email: feit1@mskcc.org
#' 
#' @examples 
#' 
#' set.seed(23420)
#' 
#' # Continuous outcome
#' dat <- simu(n=50,p=100,model="linear")
#' fit <- FLORAL(dat$xcount,dat$y,family="gaussian",progress=FALSE)
#' 
#' # Binary outcome
#' dat <- simu(n=50,p=100,model="binomial")
#' fit <- FLORAL(dat$xcount,dat$y,family="binomial",progress=FALSE)
#' 
#' # Survival outcome
#' dat <- simu(n=50,p=100,model="cox")
#' fit <- FLORAL(dat$xcount,survival::Surv(dat$t,dat$d),family="cox",progress=FALSE)
#' 
#' @import Rcpp RcppArmadillo ggplot2 RcppProgress survival glmnet dplyr grDevices utils stats
#' @importFrom survcomp concordance.index
#' @importFrom reshape melt
#' @useDynLib FLORAL
#' @export

FLORAL <- function(x,
                   y,
                   family="gaussian",
                   longitudinal=FALSE,
                   id=NULL,
                   tobs=NULL,
                   failcode=NULL,
                   length.lambda=100,
                   lambda.min.ratio=NULL,
                   mu=1,
                   ncv=5,
                   intercept=FALSE,
                   foldid=NULL,
                   step2=FALSE,
                   progress=TRUE,
                   plot=TRUE){
  
  x <- log(x+1)
  
  if (family == "gaussian"){
    
    if (longitudinal){
      stop("Longitudinal data matrix is not supported for `family=gaussian`.")
    }else{
      res <- LogRatioLasso(x,
                           y,
                           length.lambda,
                           lambda.min.ratio,
                           mu,
                           ncv,
                           intercept,
                           foldid,
                           step2,
                           progress,
                           plot)
    }
    
  }else if(family == "binomial"){
    
    if (longitudinal){
      stop("Longitudinal data matrix is not supported for `family=binomial`.")
    }else{
      res <- LogRatioLogisticLasso(x,
                                   y,
                                   length.lambda,
                                   lambda.min.ratio,
                                   mu,
                                   ncv,
                                   foldid,
                                   step2,
                                   progress,
                                   plot)
    }
    
  }else if(family == "cox"){
    
    if (longitudinal){
      
      if (is.null(id) | is.null(tobs)){
        
        stop("`id` and `tobs` must be specified if `Longitudinal` is TRUE.")
        
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
                                  length.lambda,
                                  lambda.min.ratio,
                                  mu,
                                  ncv,
                                  foldid,
                                  step2,
                                  progress,
                                  plot)
        
      }
      
    }else{
      
      res <- LogRatioCoxLasso(x,
                              y,
                              length.lambda,
                              lambda.min.ratio,
                              mu,
                              ncv,
                              foldid,
                              step2,
                              progress,
                              plot)
      
    }
  }else if (family == "finegray"){
    
    if (longitudinal){
      
      if (is.null(id) | is.null(tobs)){
        
        stop("`id` and `tobs` must be specified if `Longitudinal` is TRUE.")
        
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
          df_FG <- finegray(Surv(tstart, tstop, dd) ~ ., id=id, data=xy, etype=failcode,timefix = FALSE)
        }else{
          df_FG <- finegray(Surv(tstart, tstop, dd) ~ ., id=id, data=xy, etype=failcode,timefix = FALSE)
        }
        
        newx <- as.matrix(df_FG[,colnames(x)])
        newy <- Surv(df_FG$fgstart,df_FG$fgstop,df_FG$fgstatus)
        weight <- df_FG$fgwt
        newid <- df_FG$id
        
        res <- LogRatioFGLasso(newx,
                               newy,
                               newid,
                               weight,
                               length.lambda,
                               lambda.min.ratio,
                               mu,
                               ncv,
                               foldid,
                               step2,
                               progress,
                               plot)
        
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
                             length.lambda,
                             lambda.min.ratio,
                             mu,
                             ncv,
                             foldid,
                             step2,
                             progress,
                             plot)
      
    }
    
  }
  
}