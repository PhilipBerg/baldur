#' Spiked-in data set of reversibly oxidized cysteins
#'
#' A dataset containing quantification of reversibly oxidized cysteines using Progenesis.
#' True positives cysteines spiked-in from yeast at two different concentrations
#' and true negatives from \emph{Chlamydomonas reinhardtii} with the same
#' concentration in all samples. For details see
#' \insertCite{berg2019evaluation;textual}{baldur} and if you use this dataset
#' please cite the same paper.
#' @format A data frame with 2235 rows and 7 variables:
#' \describe{
#'   \item{identifier}{id column for features, true positives contains YEAST and
#'   true negatives contains Cre}
#'   \item{ng50_* }{Biological replicates with true positives spiked-in from 50
#'   ng yeast cells}
#'   \item{ng100_* }{Biological replicates with true positives spiked-in from 100
#'   ng yeast cells}
#' }
#' @source \url{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2619-6}
#' @references
#' \insertAllCited{}
"yeast"
