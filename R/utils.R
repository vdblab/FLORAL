
#' Clean covariate Columns
#'.This converts categorical data to 1-hot encoded columns, and retains numeric column
#' @param df Data.frame of covariates for this dataset
#' @param cols columns to clean; columns not mentioned here are dropped
#' @param drop_first default TRUE; drop first level of encoding for 1-hot
#'
#' @return dataframe of recoded covariates 
#'
#' @examples
#' df = data.frame(x=c("A", "B","C", "A", "B"))
#' clean_covariate_columns(df,cols="x")
#' clean_covariate_columns(df,cols="x", drop_first=FALSE)
#'
clean_covariate_columns <- function(df, cols, drop_first = TRUE){
  result <- data.frame(placeholder = character(nrow(df)))
  for (col in cols) {
    # keep numeric asis
    if (is.numeric(df[, col])){
      result = cbind(result, df[, col])
      colnames(result)[-1] <- col
    } else {
      vals <- unique(df[, col])
      valnames = paste0(col, "_", vals)
      if (drop_first) {
        vals <- vals[-1]
        valnames = valnames[-1]
      }
      res <- data.frame(ifelse(vals[1]  == df[, col], 1, 0))
      if (length(vals) > 1){
        for (val_i in 2:length(vals)){
          res = cbind(res, data.frame(ifelse(vals[val_i]  == df[, col], 1, 0)))
        }
      }
      colnames(res) = valnames 
      result <- cbind(result, res)
    }
  }
  #drop placeholder
  result <- result[, !colnames(result) == "placeholder"]
  # deal with single covariates where results becomes a vector
  if (!is.data.frame(result)){
    result = data.frame(result)
    colnames(result) <- cols
  }
  return(result)
}

#' Create data input list from phyloseq object
#'
#' @param phy Phyloseq object
#' @param y Outcome column of interest from phy's sample_data
#' @param covariates Covariate column names from phy's sample_data
#'
#' @return list
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' # add a covariate
#' sample_data(GlobalPatterns)$test <- rep(c(1, 0), nsamples(GlobalPatterns)/2)
#' GlobalPatterns <- tax_glom(GlobalPatterns, "Family")
#' dat <- phy_to_floral_data(GlobalPatterns, y = "test", covariates = c("SampleType"))
#' res <- FLORAL(x = dat$xcount, y=dat$y, ncov =dat$ncov, family = "Gaussian")

phy_to_floral_data<- function(phy, y=NULL, covariates=NULL){
  
  xcount = phyloseq::otu_table(phy) 
  if (nrow(xcount) != phyloseq::nsamples(phy)){
    # support both phyloseq objects with taxa as rows or columns
    xcount = t(xcount)
  }
  if(any(rownames(xcount) != phyloseq::sample_names(phy))){
    stop("malformed phyloseq object; columns of otu_table do not match sample IDs")
  }
  ncov = 0
  sampdat = phyloseq::sample_data(phy) %>% data.frame() 
  yres = sampdat %>% pull(all_of(y))
  if (!missing(covariates)){
    cov_df <- sampdat %>% select(all_of(covariates))
    cov_df_clean <- clean_covariate_columns(df=cov_df, cols = covariates)
    ncov = ncol(cov_df_clean)
    # as.matrix here speeds things significantly
    xcount = cbind(as.matrix(cov_df_clean), xcount )
  }
  return(list("xcount" = xcount, "ncov"=ncov, "y"=yres))
  
}

