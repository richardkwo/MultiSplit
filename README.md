# MultiSplit
**Hypothesis testing with multiple data splits and exchangeable p-values**

Statistical tests are sometimes constructed with data splitting. When such tests are applied to data, the result can depend on the way the data is split, which is typically random. Therefore, on a dataset, the result     of a test is random and not replicable. Further, such tests typically have low power because the full sample is not utlized. 

R package `MultiSplit` properly aggregates the results from multiple data splits and reports the p-value of the aggregated statistic. The constructed test has level that asymptotically approaches the nominal level. Typically, by aggregating results from a sufficiently large number of data splits, the test becomes replicable and much more powerful. This package implements a generic method that handles any test that is constructed with "extra randomness", including random data splitting, resampling, imputation, etc.

## Installation

The package can be installed from GitHub.

``` r
# install.packages("devtools")
# devtools::install_github("cran/kedd")
devtools::install_github("richardkwo/MultiSplit")
```

In case of problem, first make sure dependency [kedd](https://cran.r-project.org/package=kedd) is properly installed. It is removed from CRAN but can be found from CRAN [archive](https://cran.r-project.org/src/contrib/Archive/kedd/) or [GitHub](https://github.com/cran/kedd).

```R
devtools::install_github("cran/kedd")
```

For a quick start, check out the vignette in folder `vignette/` and the documentation for main function `test.multisplit`.

## Parallelization

The package supports parallelized subsampling with package [future](https://cran.r-project.org/package=future). The mode of parallelization is controlled by command `plan()` (e.g., `plan(multisession)`). See also [here](https://cran.r-project.org/web/packages/future/vignettes/future-1-overview.html). 

## Reproduction scripts

Folder `scripts/` contains code for reproducing numerical experiments. 

## Changelog

* v0.1.1
  * Fixed malfunction of `kedd` due to new `if` [behavior in R 4.2+](https://stackoverflow.com/questions/72848442/r-warning-lengthx-2-1-in-coercion-to-logical1)
* v0.1.0
  * Intial release


## Reference

Guo, F. Richard, and Rajen D. Shah. "Rank-transformed subsampling: inference for multiple data splitting and exchangeable p-values." [*arXiv:2301.02739*](https://arxiv.org/abs/2301.02739) (2023).