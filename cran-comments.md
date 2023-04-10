# For R version 4.2.3
Currently, Baldur will only be available for R versions 4.2.0 and up.

`devtools::check(remote = TRUE, manual = TRUE)`
## R CMD check results -- baldur 0.0.1 --
Duration: 1h 16m 56.1s

❯ checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Philip Berg <pb1015@msstate.edu>'
  
  New submission

❯ checking installed package size ... NOTE
    installed size is  5.4Mb
    sub-directories of 1Mb or more:
      data   1.2Mb
      libs   3.5Mb

❯ checking for future file timestamps ... NOTE
  unable to verify current time

❯ checking dependencies in R code ... NOTE
  Namespaces in Imports field not imported from:
    'RcppParallel' 'rstantools'
    All declared Imports should be used.

❯ checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

0 errors ✔ | 0 warnings ✔ | 5 notes ✖


I believe my notes are mainly from my dependencies on `rstan` and are created by `rstantools`.
As such, I hope that they are correct. The only other note is for the data.
I have two datasets for two different vignettes, they've been "optimally" compressed with `tools::resaveRdaFiles("data/")`.

In addition I have ran:
`devtools::check_win_devel()`
`devtools::check_win_release()`
which only returned one note:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Philip Berg <pb1015@msstate.edu>'

New submission

Possibly misspelled words in DESCRIPTION:
  Proteomics (2:54)
  proteomics (8:38)
    
Which are both correctly spelled.

* This is a new release.
