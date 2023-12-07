# Export the "vcfreader" C++ class by explicitly requesting vcfreader be
# exported via roxygen2's export tag.
#' @export vcfreader
#' @export vcfwriter
loadModule("vcfreader", TRUE)
loadModule("vcfwriter", TRUE)
