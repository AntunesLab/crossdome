##' Class "xrResult"
##' Description
##'
##'
##' @name xrResult-class
##' @aliases xrResult-class
##'
##' @docType class
##' @slot result Description
##' @slot allele Allele
##' @slot position_weight numeric
##' @exportClass xrResult
##' @keywords classes

# http://adv-r.had.co.nz/S4.html

setClass("xrResult",
         representation =
           representation(
             result = "data.frame",
             allele = "character",
             position_weight = "numeric",
             timestamp   = "character"
         )
)

##' Class "xrBackground"
##' Description
##'
##' @name xrBackground-class
##' @aliases xrBackground-class
##'
##' @docType class
##' @slot result Description
##' @slot allele Allele
##' @slot position_weight numeric
##' @exportClass xrBackground
##' @keywords classes

setClass("xrBackground",
         representation =
           representation(
             allele = "character",
             peptides   = "character"
           )
)
