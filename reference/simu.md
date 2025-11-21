# Simulate data following log-ratio model

Simulate a dataset from log-ratio model.

## Usage

``` r
simu(
  n = 100,
  p = 200,
  model = "linear",
  weak = 4,
  strong = 6,
  weaksize = 0.125,
  strongsize = 0.25,
  pct.sparsity = 0.5,
  rho = 0,
  timedep_slope = NULL,
  timedep_cor = NULL,
  geetype = "gaussian",
  m = 4,
  corstr = "exchangeable",
  sdvec = NULL,
  rhogee = 0.8,
  geeslope = 2.5,
  longitudinal_stability = TRUE,
  ncov = 0,
  betacov = 0,
  intercept = FALSE
)
```

## Arguments

- n:

  An integer of sample size

- p:

  An integer of number of features (taxa).

- model:

  Type of models associated with outcome variable, can be "linear",
  "binomial", "cox", "finegray", "gee" (scalar outcome with
  time-dependent features), or "timedep" (survival endpoint with
  time-dependent features).

- weak:

  Number of features with `weak` effect size.

- strong:

  Number of features with `strong` effect size.

- weaksize:

  Actual effect size for `weak` effect size. Must be positive.

- strongsize:

  Actual effect size for `strong` effect size. Must be positive.

- pct.sparsity:

  Percentage of zero counts for each sample.

- rho:

  Parameter controlling the correlated structure between taxa. Ranges
  between 0 and 1.

- timedep_slope:

  If `model` is "timedep", this parameter specifies the slope for the
  feature trajectories. Please refer to the Simulation section of the
  manuscript for more details.

- timedep_cor:

  If `model` is "timedep", this parameter specifies the sample-wise
  correlations between longitudinal features. Please refer to the
  Simulation section of the manuscript for more details.

- geetype:

  If `model` is "gee", `geetype` is the type of GEE outcomes. Now
  support "gaussian" and "binomial".

- m:

  If `model` is "gee", `m` is the number of repeated measurements per
  subject.

- corstr:

  If `model` is "gee", `corstr` is the working correlation structure.
  Now support "independence", "exchangeable", and "AR-1".

- sdvec:

  If `model` is "gee" and `geetype` is "gaussian", `sdvec` is the vector
  of standard deviations of each outcome variable.

- rhogee:

  If `model` is "gee", `rhogee` is the correlation parameter between
  longitudinal outcomes under the selected working correlation
  structure.

- geeslope:

  If `model` is "gee", `geeslope` is the linear time effect.

- longitudinal_stability:

  If `model` is "timedep", this is a binary indicator which determines
  whether the trajectories are more stable (`TRUE`) or more volatile
  (`FALSE`).

- ncov:

  Number of covariates that are not compositional features.

- betacov:

  Coefficients corresponding to the covariates that are not
  compositional features.

- intercept:

  Boolean. If TRUE, then a random intercept will be generated in the
  model. Only works for `linear` or `binomial` models.

## Value

A list with simulated count matrix `xcount`, log1p-transformed count
matrix `x`, outcome (continuous `y`, continuous centered `y0`, binary
`y`, or survival `t`, `d`), true coefficient vector `beta`, list of
non-zero features `idx`, value of intercept `intercept` (if applicable).

## References

Fei T, Funnell T, Waters N, Raj SS et al. Enhanced Feature Selection for
Microbiome Data using FLORAL: Scalable Log-ratio Lasso Regression
bioRxiv 2023.05.02.538599.

## Author

Teng Fei. Email: feit1@mskcc.org

## Examples

``` r
set.seed(23420)
dat <- simu(n=50,p=30,model="linear")
```
