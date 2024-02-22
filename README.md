
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

![cran-badge](http://www.r-pkg.org/badges/version/baldur)
![R-CMD-check](https://github.com/PhilipBerg/baldur/actions/workflows/check-standard.yaml/badge.svg)
![Static
Badge](https://img.shields.io/badge/10.1016%2Fj.mcpro.2023.100658-green?style=palstic&label=DOI%3A&labelColor=grey&link=https%3A%2F%2Fdoi.org%2F10.1016%2Fj.mcpro.2023.100658)
<!-- badges: end -->

# Baldur

The goal of Baldur is to shine light on your proteomics data. Baldur is
a hierarchical Bayesian model that uses an empirical Bayes method to
estimate hyperparamters for the variance and measurement specific
uncertainty. It then estimates the posterior of the difference in means
between different conditions for each peptide/protein/PTM. Finally, it
integrates the posterior to estimate the probability of error.

## Installation

You can install the development version of `baldur` from this github or
the stable version from CRAN. Importantly, you first need to follow the
instructions for installing `rstan`
<https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started> and make
sure that is working. Then you can install `baldur` accordingly:

For the stable release please install from CRAN:

``` r
install.packages('baldur')
```

Or you can install the developmental version of `baldur` from this
github:

``` r
devtools::install_github('PhilipBerg/baldur', build_vignettes = T)
```

Note that Ubuntu operating systems can require `pandoc`
<https://pandoc.org/> to compile the vignettes.

For Windows, the developmental version of `rstan` is sometimes needed to
install `baldur`.

## Example

Please see the vignettes for examples
`vignette('baldur_yeast_tutorial')` and
`vignette('baldur_ups_tutorial')`.

## Reference

Berg, Philip, and George Popescu. “Baldur: Bayesian Hierarchical
Modeling for Label-Free Proteomics with Gamma Regressing Mean-Variance
Trends” Molecular & Cellular Proteomics (2023): 2023-12.
<https://doi.org/10.1016/j.mcpro.2023.100658>
