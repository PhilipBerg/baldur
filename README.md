
<!-- README.md is generated from README.Rmd. Please edit that file -->

![cran-badge](http://www.r-pkg.org/badges/version/baldur)

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
sure that is working. Then you can install `baldur` accordingly: For the
stable release please install from CRAN:

``` r
install.packages('baldur') # Once on CRAN
#> Installing package into 'C:/Users/Pberg/AppData/Local/Temp/RtmpkHUppW/temp_libpath17381ca63c17'
#> (as 'lib' is unspecified)
#> Warning: package 'baldur' is not available for this version of R
#> 
#> A version of this package for your version of R might be available elsewhere,
#> see the ideas at
#> https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages
```

Or you can install the development version of `baldur` from this github:

``` r
devtools::install_github('PhilipBerg/baldur', build_vignettes = T)
#> Downloading GitHub repo PhilipBerg/baldur@HEAD
#> cli        (3.4.1 -> 3.6.1) [CRAN]
#> rlang      (1.0.6 -> 1.1.0) [CRAN]
#> ps         (1.7.2 -> 1.7.4) [CRAN]
#> vctrs      (0.5.2 -> 0.6.1) [CRAN]
#> tibble     (3.1.8 -> 3.2.1) [CRAN]
#> pillar     (1.8.1 -> 1.9.0) [CRAN]
#> gtable     (0.3.1 -> 0.3.3) [CRAN]
#> loo        (2.5.1 -> 2.6.0) [CRAN]
#> ggplot2    (3.4.1 -> 3.4.2) [CRAN]
#> dplyr      (1.1.0 -> 1.1.1) [CRAN]
#> multidplyr (0.1.2 -> 0.1.3) [CRAN]
#> rstantools (2.3.0 -> 2.3.1) [CRAN]
#> Skipping 2 packages ahead of CRAN: StanHeaders, rstan
#> Installing 12 packages: cli, rlang, ps, vctrs, tibble, pillar, gtable, loo, ggplot2, dplyr, multidplyr, rstantools
#> Installing packages into 'C:/Users/Pberg/AppData/Local/Temp/RtmpkHUppW/temp_libpath17381ca63c17'
#> (as 'lib' is unspecified)
#> 
#>   There are binary versions available but the source versions are later:
#>         binary source needs_compilation
#> ps       1.7.3  1.7.4              TRUE
#> ggplot2  3.4.1  3.4.2             FALSE
#> 
#> package 'cli' successfully unpacked and MD5 sums checked
#> package 'rlang' successfully unpacked and MD5 sums checked
#> package 'vctrs' successfully unpacked and MD5 sums checked
#> package 'tibble' successfully unpacked and MD5 sums checked
#> package 'pillar' successfully unpacked and MD5 sums checked
#> package 'gtable' successfully unpacked and MD5 sums checked
#> package 'loo' successfully unpacked and MD5 sums checked
#> package 'dplyr' successfully unpacked and MD5 sums checked
#> package 'multidplyr' successfully unpacked and MD5 sums checked
#> package 'rstantools' successfully unpacked and MD5 sums checked
#> 
#> The downloaded binary packages are in
#>  C:\Users\Pberg\AppData\Local\Temp\RtmpMlchlD\downloaded_packages
#> installing the source packages 'ps', 'ggplot2'
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>          checking for file 'C:\Users\Pberg\AppData\Local\Temp\RtmpMlchlD\remotes46b0613abb7\PhilipBerg-baldur-d649aff/DESCRIPTION' ...  ✔  checking for file 'C:\Users\Pberg\AppData\Local\Temp\RtmpMlchlD\remotes46b0613abb7\PhilipBerg-baldur-d649aff/DESCRIPTION'
#>       ─  preparing 'baldur': (352ms)
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   ✔  checking DESCRIPTION meta-information
#> ─  cleaning src
#>       ─  installing the package to process help pages
#>      Loading required namespace: baldur
#>       ─  saving partial Rd database (1.2s)
#>          creating vignettes ...     creating vignettes ...   ✔  creating vignettes (1.4s)
#>       ─  cleaning src
#>       ─  checking for LF line-endings in source and make files and shell scripts
#>       ─  checking for empty or unneeded directories
#>       ─  building 'baldur_0.0.1.tar.gz'
#>   Warning:     Warning: file 'baldur/configure' did not have execute permissions: corrected
#>      
#> 
#> Installing package into 'C:/Users/Pberg/AppData/Local/Temp/RtmpkHUppW/temp_libpath17381ca63c17'
#> (as 'lib' is unspecified)
```

Note that Ubuntu operating systems can require `pandoc`
<https://pandoc.org/> to compile the vignettes.

## Example

Please see the vignettes for an examples
`vignette('baldur_yeast_tutorial')` and
`vignette('baldur_ups_tutorial')`.
