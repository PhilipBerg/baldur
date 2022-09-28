`devtools::check(remote = TRUE, manual = TRUE)`
## R CMD check results -- baldur 0.0.1 --
Duration: 14m 42.9s

> checking PDF version of manual ... WARNING
  LaTeX errors when creating PDF version.
  This typically indicates Rd problems.
  LaTeX errors found:
  !pdfTeX error: pdflatex (file ts1-zi4r): Font ts1-zi4r at 600 not found
   ==> Fatal error occurred, no output PDF file produced!

> checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Philip Berg <pb1015@msstate.edu>'
  
  New submission

> checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

> checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    '-mmmx' '-msse' '-msse3' '-msse4.1' '-msse4.2' '-mssse3'

> checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    'baldur-manual.tex'

0 errors v | 1 warning x | 4 notes x

The warning, I think, is due to a local issue on machine when compiling pdf documents.
Because, I have ran:
`devtools::check_win_devel()`
`rhub::check_for_cran()`
`rhub::check(platform = 'ubuntu-rchk')`
`rhub::check_with_sanitizers()`
without any problems.
I believe my notes are mainly from my dependencies on `rstan` and are created by `rstantools`.
As such, I hope that they are correct.


* This is a new release.
