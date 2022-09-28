
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Baldur

The goal of Baldur is to shine light on your proteomics data.

## Installation

Please note that `baldur` is currently only working on R versions 4.0.X.
This is due to the dependency on RStan that is not yet working on the
new toolchain associated with R versions 4.1/2.X. We hope to fix this
problem as soon as possible.

You can install the development version of `baldur` from this github.

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
if (!require("remotes")) {
  install.packages("remotes")
}
remotes::install_github("stan-dev/rstantools")
devtools::install_github('PhilipBerg/baldur', build_vignettes = T)
```

Please note that you need to follow the guidelines for setting up RStan
according to their homepage.

## Example

Please see the vignette for an example `vignette('baldur-tutorial')`.
