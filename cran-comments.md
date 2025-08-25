# For release of Baldur 0.0.4
* Added capability to add auxiliary variables to the design matrix, e.g., batches.
* Fixed bug with numeric overflow when checking that the contrast matrix is correctly setup.
* Updated links to other packages as requested
* Improved Vignettes


# For release of Baldur 0.0.3
Updated Stan models for future comparability with `rstan` as suggested by one of their developers on `baldur`'s GitHub:
https://github.com/PhilipBerg/baldur/pull/16
All checks pass with same or fewer notes reported for earlier releases.
The main difference is a note for recent update six days ago, but this update is due to unforeseen updates in the `rstan` dependency.


# For release of Baldur 0.0.2
I've patched minor bug in a function that checks that input structure is correct for several functions.
I've updated the printing of a S3 method to use Greek letters to better link the to the manuscript of `baldur` .
I've also updated one of the Stan models in `baldur` to accommodate peer-reviewer feedback.


# For R version 4.2.3
Currently, Baldur will only be available for R versions 4.2.0 and up.

`devtools::check(remote = TRUE, manual = TRUE)`
## R CMD check results -- baldur 0.0.1.1 --
Duration: 48m 23.2s

❯ checking installed package size ... NOTE
    installed size is  5.4Mb
    sub-directories of 1Mb or more:
      data   1.2Mb
      libs   3.5Mb

❯ checking dependencies in R code ... NOTE
  Namespaces in Imports field not imported from:
    'RcppParallel' 'rstantools'
    All declared Imports should be used.

❯ checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

0 errors ✔ | 0 warnings ✔ | 3 notes ✖


I believe my notes are mainly from my dependencies on `rstan` and are created by `rstantools`.
As such, I hope that they are correct. The only other note is for the data.
I have two datasets for two different vignettes, they've been "optimally" compressed with `tools::resaveRdaFiles("data/")`.

In addition I have ran:
`devtools::check_win_devel()`
`devtools::check_win_release()`
which only returned one note:

* checking CRAN incoming feasibility ... [10s] NOTE
Maintainer: 'Philip Berg <pb1015@msstate.edu>'

Possibly misspelled words in DESCRIPTION:
  Popescu (18:18)
    
Which is correctly spelled.

Finally, my GitHub workflows for R-CMD-check passed on Ubuntu and MacOS.
