
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Baldur

The goal of Baldur is to shine light on your proteomics data.

## Installation

You can install the development version of `baldur` from this github or
the stable version from CRAN. Importantly, you first need to follow the
instructions for installing `rstan`
<https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started>. Then you
can install `baldur` as accordingly:

You can install the development version of `baldur` from this GitHub:

```{r}
devtools::install_github('PhilipBerg/baldur', build_vignettes = T)
# or
install.package('baldur') # Once on CRAN
```

## Example

Please see the vignettes for an examples
`vignette('baldur_yeast_tutorial')` and
`vignette('baldur_ups_tutorial')`.
