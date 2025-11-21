# Fit Log-ratio lasso regression for compositional covariates

Conduct log-ratio lasso regression for continuous, binary and survival
outcomes.

## Usage

``` r
FLORAL(
  x,
  y,
  ncov = 0,
  family = "gaussian",
  longitudinal = FALSE,
  id = NULL,
  tobs = NULL,
  failcode = NULL,
  corstr = "exchangeable",
  scalefix = FALSE,
  scalevalue = 1,
  pseudo = 1,
  length.lambda = 100,
  lambda.min.ratio = NULL,
  ncov.lambda.weight = 0,
  a = 1,
  mu = 1,
  pfilter = 0,
  maxiter = 100,
  ncv = 5,
  ncore = 1,
  intercept = FALSE,
  foldid = NULL,
  step2 = TRUE,
  progress = TRUE,
  plot = TRUE
)
```

## Arguments

- x:

  Feature matrix, where rows specify subjects and columns specify
  features. The first `ncov` columns should be patient characteristics
  and the rest columns are microbiome absolute counts corresponding to
  various taxa. If `x` contains longitudinal data, the rows must be
  sorted in the same order of the subject IDs used in `y`.

- y:

  Outcome. For a continuous or binary outcome, `y` is a vector. For
  survival outcome, `y` is a `Surv` object.

- ncov:

  An integer indicating the number of first `ncov` columns in `x` that
  will not be subject to the zero-sum constraint.

- family:

  Available options are `gaussian`, `binomial`, `cox`, `finegray`.

- longitudinal:

  `TRUE` or `FALSE`, indicating whether longitudinal data matrix is
  specified for input `x`. (`Longitudinal=TRUE` and `family="cox"` or
  `"finegray"` will fit a time-dependent covariate model.
  `Longitudinal=TRUE` and `family="gaussian"` or `"binomial"` will fit a
  GEE model.)

- id:

  If `longitudinal` is `TRUE`, `id` specifies subject IDs corresponding
  to the rows of input `x`.

- tobs:

  If `longitudinal` is `TRUE`, `tobs` specifies time points
  corresponding to the rows of input `x`.

- failcode:

  If `family = finegray`, `failcode` specifies the failure type of
  interest. This must be a positive integer.

- corstr:

  If a GEE model is specified, then `corstr` is the corresponding
  working correlation structure. Options are `independence`,
  `exchangeable`, `AR-1` and `unstructured`.

- scalefix:

  `TRUE` or `FALSE`, indicating whether the scale parameter is estimated
  or fixed if a GEE model is specified.

- scalevalue:

  Specify the scale parameter if `scalefix=TRUE`.

- pseudo:

  Pseudo count to be added to `x` before taking log-transformation. If
  unspecified, then the log-transformation will not be performed.

- length.lambda:

  Number of penalty parameters used in the path

- lambda.min.ratio:

  Ratio between the minimum and maximum choice of lambda. Default is
  `NULL`, where the ratio is chosen as 1e-2.

- ncov.lambda.weight:

  Weight of the penalty lambda applied to the first `ncov` covariates.
  Default is 0 such that the first `ncov` covariates are not penalized.

- a:

  A scalar between 0 and 1: `a` is the weight for lasso penalty while
  `1-a` is the weight for ridge penalty.

- mu:

  Value of penalty for the augmented Lagrangian

- pfilter:

  A pre-specified threshold to force coefficients with absolute values
  less than pfilter times the maximum value of absolute coefficient as
  zeros in the GEE model. Default is zero, such that all coefficients
  will be reported.

- maxiter:

  Number of iterations needed for the outer loop of the augmented
  Lagrangian algorithm.

- ncv:

  Folds of cross-validation. Use `NULL` if cross-validation is not
  wanted.

- ncore:

  Number of cores for parallel computing for cross-validation. Default
  is 1.

- intercept:

  `TRUE` or `FALSE`, indicating whether an intercept should be
  estimated.

- foldid:

  A vector of fold indicator. Default is `NULL`.

- step2:

  `TRUE` or `FALSE`, indicating whether a second-stage feature selection
  for specific ratios should be performed for the features selected by
  the main lasso algorithm. Will only be performed if cross validation
  is enabled.

- progress:

  `TRUE` or `FALSE`, indicating whether printing progress bar as the
  algorithm runs.

- plot:

  `TRUE` or `FALSE`, indicating whether returning plots of model
  fitting.

## Value

A list with path-specific estimates (beta), path (lambda), and others.
Details can be found in `README.md`.

## References

Fei T, Funnell T, Waters N, Raj SS et al. Enhanced Feature Selection for
Microbiome Data using FLORAL: Scalable Log-ratio Lasso Regression
bioRxiv 2023.05.02.538599.

## Author

Teng Fei. Email: feit1@mskcc.org

## Examples

``` r
set.seed(23420)

# Continuous outcome
dat <- simu(n=50,p=30,model="linear")
fit <- FLORAL(dat$xcount,dat$y,family="gaussian",ncv=2,progress=FALSE,step2=TRUE)

# Binary outcome
# dat <- simu(n=50,p=30,model="binomial")
# fit <- FLORAL(dat$xcount,dat$y,family="binomial",progress=FALSE,step2=TRUE)

# Survival outcome
# dat <- simu(n=50,p=30,model="cox")
# fit <- FLORAL(dat$xcount,survival::Surv(dat$t,dat$d),family="cox",progress=FALSE,step2=TRUE)

# Competing risks outcome
# dat <- simu(n=50,p=30,model="finegray")
# fit <- FLORAL(dat$xcount,survival::Surv(dat$t,dat$d,type="mstate"),failcode=1,
#               family="finegray",progress=FALSE,step2=FALSE)

# Longitudinal continuous outcome
#  dat <- simu(n=50,p=30,model="gee",geetype="gaussian",m=3,corstr="exchangeable",sdvec=rep(1,3))
#  fit <- FLORAL(x=cbind(dat$tvec, dat$xcount),y=dat$y,id=dat$id,family="gaussian",
#                ncov=1,longitudinal = TRUE,corstr = "exchangeable",lambda.min.ratio=1e-3,
#                progress=FALSE,step2=FALSE)
```
