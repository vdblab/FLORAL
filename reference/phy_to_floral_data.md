# Create data input list from phyloseq object

Create data input list from phyloseq object

## Usage

``` r
phy_to_floral_data(phy, y = NULL, covariates = NULL)
```

## Arguments

- phy:

  Phyloseq object

- y:

  Outcome column of interest from phy's sample_data

- covariates:

  Covariate column names from phy's sample_data

## Value

list

## Examples

``` r
library(phyloseq)
data(GlobalPatterns)
# add a covariate
sample_data(GlobalPatterns)$test <- rep(c(1, 0), nsamples(GlobalPatterns)/2)
# GlobalPatterns <- tax_glom(GlobalPatterns, "Phylum")
dat <- phy_to_floral_data(GlobalPatterns, y = "test", covariates = c("SampleType"))
# res <- FLORAL(x = dat$xcount, y=dat$y, ncov=dat$ncov, family = "binomial", ncv=NULL)
```
