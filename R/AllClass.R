#' Class "xrBackground"
#' @description An object holding data from \code{\link{hla_database}}
#'
#' @name xrBackground-class
#' @docType class
#'
#' @slot allele Description
#' @slot peptides Description
#' @slot stats A slot containing statistics
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
#' @slot query A peptide target. Only 9-mers are supported.
#' @slot result Crossdome ranking data.frame
#' @slot allele Input a MHC Class I allele
#' @slot expression Expression data slot
#' @slot analysis Sequence and immunogenicity slot
#' @slot position_weight A numeric vector derived from TCR hotspots
#' @slot timestamp Execution timestamp
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
