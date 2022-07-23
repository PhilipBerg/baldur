
# Baldur
The goal of Baldur is to shine light on your proteomics data.

## Installation
Please note that `baldur` is currently only working on R versions 4.0.X.
This is due to the dependency on RStan that is not yet working on the new toolchain associated with R versions 4.1.X.
We hope to fix this problem as soon as possible.

You can install the development version of `baldur` as follows:
``` r
require(devtools)
devtools::install_github('PhilipBerg/baldur', build_vignettes = T)
```

Please note that you need to follow the guidelines for setting up RStan according to their homepage.

## Example

Please see the vignette for an example `vignette('baldur-tutorial')`.
