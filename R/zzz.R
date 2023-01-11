##' @importFrom utils packageDescription
# Should I mentioned third-party software?

.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields = "Version")

  msg <- paste0(pkgname, " v", pkgVersion, ". For help, please check the package vignettes!:) \n")
  citation <- paste0("If you use ", pkgname, " in published research, please cite: <coming soon>\n")

  packageStartupMessage(paste0(msg, citation))
}
