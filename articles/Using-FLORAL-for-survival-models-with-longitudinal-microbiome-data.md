# Using FLORAL for survival models with longitudinal microbiome data

``` r
library(FLORAL)
library(dplyr)
library(patchwork)
library(survival)
set.seed(8192024)
```

In this vignette, we illustrate how to apply `FLORAL` to fit a Cox model
with longitudinal microbiome data. Due to limited availability of public
data sets with survival information, we use simulated data for
illustrative purposes.

## Data simulation

We will use the built-in simulation function
[`simu()`](https://vdblab.github.io/FLORAL/reference/simu.md) to
generate longitudinal compositional features and the corresponding
time-to-event. The underlying methodology used for the simulation is
based on a piece-wise exponential distribution as described by [Hendry
2014](https://doi.org/10.1002/sim.5945).

By default, the first 10 features out of the 500 features simulated
below are associated with the time-to-event.

``` r

simdat <- simu(n=200, # sample size
               p=500, # number of features
               model="timedep",
               pct.sparsity = 0.8, # proportion of zeros
               rho=0, # feature-wise correlation
               longitudinal_stability = TRUE # choose to simulate longitudinal features with stable trajectories
)
```

With the simulated data, the log-ratio lasso Cox model with
time-dependent features can be fitted by running the following function.
Here we provide a detailed description on each arguments:

- First of all, please use `longitudinal = TRUE` such that the algorithm
  would use the appropriate method to handle longitudinal data.
- The feature matrix input `x` should be the count matrix where rows
  specify samples and columns specify features.
- The vector of IDs of subjects/patients corresponding to the rows of
  `x` should be input as `id`.
- The vector of sample collection times corresponding to the rows of `x`
  should be input as `tobs`.
- The `Surv` object (`Surv(time,status)`) of **unique patients** should
  be input as `y`. Please note that the survival data should be sorted
  with respect to the IDs specified in `id`.

``` r

fit <- FLORAL(x=simdat$xcount,
              y=Surv(simdat$data_unique$t,simdat$data_unique$d),
              family="cox",
              longitudinal = TRUE,
              id = simdat$data$id,
              tobs = simdat$data$t0,
              progress=FALSE,
              plot=TRUE)

fit$selected
#> $min
#> [1] "taxa1"   "taxa2"   "taxa27"  "taxa366" "taxa38"  "taxa5"   "taxa6"  
#> [8] "taxa8"   "taxa9"  
#> 
#> $`1se`
#> [1] "taxa1" "taxa5" "taxa6" "taxa8" "taxa9"
#> 
#> $min.2stage
#> [1] "taxa2"   "taxa366" "taxa38"  "taxa5"   "taxa6"   "taxa8"   "taxa9"  
#> 
#> $`1se.2stage`
#> [1] "taxa1" "taxa5" "taxa6" "taxa8" "taxa9"
```

The list of selected features is saved in `fit$selected` as shown above.

To appropriately prepare the data in practice, we have the following
recommendations:

- Start with patient metadata which includes survival data (time and
  status), sorting the metadata by patient IDs. Extract time and status
  variables for the `Surv` object for input as `y`.
- Curate the microbiome feature data matrix, sorted by patient IDs and
  time of sample collection. Save the patient ID and time of sample
  collection vectors for `id` and `tobs`. Save the feature table for
  input as `x`.
