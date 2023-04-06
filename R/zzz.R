.onLoad <- function(libname, pkgname) {
  lapply(seq_along(stanmodels), \(x) assign(names(stanmodels)[x], stanmodels[[x]], envir = asNamespace('baldur')))
  Rdpack::Rdpack_bibstyles(package = pkgname, authors = "LongNames")
}
