#' Class "xrBackground"
#' Description
#'
#' @name xrBackground-class
#' @docType class
#'
#' @slot peptides Description
#' @slot allele Description
#'
#' @exportClass xrBackground
#' @keywords classes

# http://adv-r.had.co.nz/S4.html

setClass("xrBackground",
         representation =
           representation(
             allele = "character",
             peptides   = "character"
           )
)

#' Class "xrResult"
#' Description
#'
#'
#' @name xrResult-class
#' @docType class
#'
#' @slot result Description
#' @slot allele Allele
#' @slot position_weight numeric
#'
#' @exportClass xrResult
#' @keywords classes

setClass("xrResult",
         representation =
           representation(
             result = "data.frame",
             allele = "character",
             position_weight = "numeric",
             timestamp   = "character"
         )
)
