# FLORAL (development version)

* Adding instructions to install the package via `bioconda`.

* Address new CRAN requirements regarding Internet access.

# FLORAL 0.4.0

* Adding a function `phy_to_floral_data` which helps convert a `phyloseq` object to be compatible with `FLORAL`.

* Adding a corresponding vignette for `phy_to_floral_data`.

* Introducing the new GEE model for longitudinal continuous and binary outcomes. More documentations to follow in the next development version.

# FLORAL 0.3.0

* Improves stability when fitting Fine-Gray model with longitudinal covariates.

* Enables parallel computation for cross-validation.

* Fixes several bugs as reported in Issues. 

* Adding options to use user-specified pseudo counts

* Adding options to use user-specified number of maximum number of iterations

* Adding a simulation scenario for survival regression models with longitudinal features

* (BETA version of) the new GEE method

# FLORAL 0.2.0

* Including more examples in document compared to CRAN version (0.1.0.9000)

* Enables elastic net models. Users can specify the weight of lasso penalty using argument `a`. (0.1.0.9001)

* Allows adding non-compositional covariates which are not constrained by the zero-sum constraint. (0.1.0.9001)

* Adds a function `mcv.FLORAL()` to perform multiple runs of k-fold cross-validation to summarize selection probabilities for features. (0.1.0.9001)

* Adds a function `a.FLORAL()` to compare different choices of elastic net weight `a` for a fixed cross-validation setting. (0.1.0.9001)

# FLORAL 0.1.0

* Initial release.