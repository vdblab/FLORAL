#' Fit Log-ratio lasso regression for compositional covariates 
#'
#' @description Conduct log-ratio lasso regression for continuous, binary and survival outcomes. 
#' @param x Count data matrix, where rows specify subjects and columns specify features. If \code{x} contains longitudinal data, the rows must be sorted in the same order of the subject IDs used in \code{y}.
#' @param y Outcome. For a continuous or binary outcome, \code{y} is a vector. For survival outcome, \code{y} is a \code{Surv} object.
#' @param family Available options are \code{gaussian}, \code{binomial}, \code{cox}, \code{finegray}.
#' @param longitudinal \code{TRUE} or \code{FALSE}, indicating whether longitudinal data matrix is specified for input \code{x}. (Still under development. Please use with caution)
#' @param id If \code{longitudinal} is \code{TRUE}, \code{id} specifies subject IDs corresponding to the rows of input \code{x}.
#' @param tobs If \code{longitudinal} is \code{TRUE}, \code{tobs} specifies time points corresponding to the rows of input \code{x}.
#' @param failcode If \code{family = finegray}, \code{failcode} specifies the failure type of interest. This must be a positive integer.
#' @param length.lambda Number of penalty parameters used in the path
#' @param lambda.min.ratio Ratio between the minimum and maximum choice of lambda. Default is \code{NULL}, where the ratio is chosen as 1e-2.
#' @param mu Value of penalty for the augmented Lagrangian
#' @param ncv Number of cross-validation runs. Use \code{NULL} if cross-validation is not wanted.
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
#' fit <- FLORAL(dat$xcount,dat$y,family="gaussian",progress=FALSE,step2=TRUE)
#' 
#' @import Rcpp ggplot2 survival glmnet dplyr
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
                   step2=TRUE,
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
        
        stop("`id` and `tobs` must be specified if `longitudinal` is TRUE.")
        
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
        
        stop("`id` and `tobs` must be specified if `longitudinal` is TRUE.")
        
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
  
  res$loss <- NULL
  res$mse <- NULL
  res$best.beta$min <- res$best.beta$min.mse
  res$best.beta$`1se` <- res$best.beta$add.1se
  res$best.beta$min.mse <- NULL
  res$best.beta$add.1se <- NULL
  
  if (step2){
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
      res$step2.ratios$min.idx <- res$step2.feature.min[,!is.na(colSums(res$step2.feature.min))]
      
      res$step2.tables$`min` <- summary(res$step2fit.min)$coefficients
      res$step2.tables$`min` <- res$step2.tables$`min`[rownames(res$step2.tables$`min`) != "(Intercept)", ]
      rownames(res$step2.tables$`min`) <- res$step2.ratios$`min`
    }
    if (length(res$selected.feature$`1se.2stage`)>0){
      namemat <- matrix(names(res$best.beta$`1se`)[res$step2.feature.1se],nrow=2)
      res$step2.ratios$`1se` <- as.vector(na.omit(apply(namemat,2,function(x) ifelse(sum(is.na(x))==0,paste(x,collapse ="/"),NA))))
      res$step2.ratios$`1se.idx` <- res$step2.feature.1se[,!is.na(colSums(res$step2.feature.1se))]
      
      res$step2.tables$`1se` <- summary(res$step2fit.1se)$coefficients
      res$step2.tables$`1se` <- res$step2.tables$`1se`[rownames(res$step2.tables$`1se`) != "(Intercept)", ]
      rownames(res$step2.tables$`1se`) <- res$step2.ratios$`1se`
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