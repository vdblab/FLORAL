
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FLORAL: Log-ratio Lasso Regression

<!-- badges: start -->

[![R-CMD-check](https://github.com/vdblab/FLORAL/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vdblab/FLORAL/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`FLORAL` is an open-source computational tool to perform log-ratio lasso
regression modeling and microbial feature selection for continuous,
binary, time-to-event, and competing risk outcomes. The proposed method
adapts the augmented Lagrangian algorithm for a zero-sum constraint
optimization problem while enabling a two-stage screening process for
extended false-positive control.

## Installation

You can install the development version of `FLORAL` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("vdblab/FLORAL")
```

## Example

Here is a toy example for linear regression with 10-fold
cross-validation for a simulated data with 50 samples and 100
compositional features.

The data simulation procedure is described in the preprint.

``` r

set.seed(23420)
library(FLORAL)

dat <- simu(n=50,p=100,model="linear")
fit <- FLORAL(dat$xcount,dat$y,family="gaussian",ncv=10,progress=FALSE)
```

To view plots of cross-validated prediction error and parameter
coefficients, use `fit$pmse` or `fit$pcoef`:

<img src="man/figures/README-plot-1.png" width="100%" /><img src="man/figures/README-plot-2.png" width="100%" />

For binary and survival outcomes, please specify `family="binomial"` and
`family="cox"` accordingly.

``` r

library(FLORAL)
library(survival)
#> Warning: package 'survival' was built under R version 4.1.3

dat.bin <- simu(n=50,p=100,model="binomial")
fit.bin <- FLORAL(dat.bin$xcount,dat.bin$y,family="binomial",ncv=10,progress=FALSE)

dat.cox <- simu(n=50,p=100,model="cox")
fit.cox <- FLORAL(dat.cox$xcount,Surv(dat.cox$t,dat.cox$d),family="cox",ncv=10,progress=FALSE)
```
