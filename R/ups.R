#' Spiked-in data set of peptides
#'
#' A dataset containing quantification of peptides using Progenesis. True
#' positives peptides spiked-in from the Universal Proteomics Standard Set 1
#' (UPS1) at three different concentrations and true negatives from
#' \emph{Chlamydomonas reinhardtii} with the same concentration in all samples.
#' You can find true positives with `stringr::str_detect(ups$identifier,
#' 'UPS')`. For details see \insertCite{berg2019evaluation;textual}{baldur} and
#' if you use this dataset please cite the same paper.
#' @format A data frame with 10599 rows and 13 variables:
#'  \describe{
#'   \item{identifier}{id column for features, true positives contains UPS and
#'   true negatives contains Cre}
#'    \item{fmol25_1}{First technical replicate
#'   with true positives spiked-in from 25 fmol UPS1 peptides}
#'    \item{fmol25_2}{Second technical replicate
#'    with true positives spiked-in from 25 fmol
#'   UPS1 peptides}
#'    \item{fmol25_3}{Third technical replicate with true
#'   positives spiked-in from 25 fmol UPS1 peptides}
#'    \item{fmol25_4}{Fourth
#'   technical replicate with true positives spiked-in from 25 fmol UPS1
#'   peptides}
#'   \item{fmol50_1}{First technical replicate with true positives
#'   spiked-in from 50 fmol UPS1 peptides}
#'    \item{fmol50_2}{Second technical
#'   replicate with true positives spiked-in from 50 fmol UPS1 peptides}
#'   \item{fmol50_3}{Third technical replicate with true positives spiked-in
#'   from 50 fmol UPS1 peptides}
#'   \item{fmol50_4}{Fourth technical replicate
#'   with true positives spiked-in from 50 fmol UPS1 peptides}
#'   \item{fmol100_1}{First technical replicate with true positives
#'   spiked-in from 100 fmol
#'   UPS1 peptides}
#'   \item{fmol100_2}{Second technical replicate with true
#'   positives spiked-in from 100 fmol UPS1 peptides}
#'   \item{fmol100_3}{Third
#'   technical replicate with true positives spiked-in from 100 fmol UPS1
#'   peptides}
#'   \item{fmol100_4}{Fourth technical replicate with true positives
#'   spiked-in from 100 fmol UPS1 peptides}
#'   }
#' @source
#'   \url{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2619-6}
#'
#' @references
#'   \insertAllCited{}
"ups"
