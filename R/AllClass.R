#' Class "xrBackground"
#' @description An object holding data from \code{\link{hla_database}}
#'
#' @name xrBackground-class
#' @docType class
#'
#' @slot peptides Description
#' @slot allele Description
#'
#' @exportClass xrBackground
#' @keywords classes

setClass("xrBackground",
         slots = c(
             allele = "character",
             peptides   = "character",
             stats = "list"
           ),
         validity = function(object) {
           ninemer <- all(sapply(object@peptides, nchar) == 9)
           if(!ninemer) {
             print("Please, Crossdome only supports 9-mer peptides.")
           }
         }
)

#' Class "xrResult"
#' @description The Crossdome main object
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
         slots = c(
             query = "character",
             result = "data.frame",
             allele = "character",
             expression = "list",
             analysis = "list",
             position_weight = "numeric",
             timestamp   = "character"
         )
)
