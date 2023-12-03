## usethis namespace: start
#' Imports
#' @useDynLib vcfppR, .registration = TRUE
#' @importFrom Rcpp sourceCpp loadModule
## usethis namespace: end
"_PACKAGE"
Rcpp::loadModule("vcfreader", TRUE)
