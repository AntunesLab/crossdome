#' .internal_checking_epitope
#' @description Returns MDS components for each amine acid
#'
#' @import Biostrings
#' @importFrom Biostrings AA_STANDARD
#' @noRd

.internal_checking_epitope <- function(epitope) {
  epitope <- base::strsplit(epitope, split = "")[[1]]
  if(!all(epitope %in% Biostrings::AA_STANDARD)) {
    quit("Please, check your input sequence.")
  }
  return(epitope)
}

#' .internal_epitope_to_matrix
#' @description Returns MDS components for each amine acid
#'
#' @import Biostrings
#' @importFrom Biostrings AA_STANDARD
#' @importFrom stats cor
#' @noRd

.internal_epitope_to_matrix <- function(epitope) {
  if(length(epitope) <= 1) {
    epitope <- .internal_checking_epitope(epitope)
  }
  epitope_to_matrix <- MDS_COMPONENTS[, epitope]
  return(epitope_to_matrix)
}

#' .internal_matrix_metrics
#' @description Returns Frobenius and Braun et al score
#' @noRd

.internal_matrix_metrics <- function(
    query_components, subject_components, position_weight = NULL) {

  product_components <- (query_components - subject_components)**2
  relatedness_score <- sqrt(apply(product_components, 2, sum))

  if(!is.null(position_weight)) {
    relatedness_score <- relatedness_score * sqrt(position_weight)
  }

  relatedness_score <- sum(relatedness_score) / ncol(query_components)

  return(list(relatedness_score = relatedness_score))
}

#' head
#' @description Print out crossdome result object
#' @noRd
head.xr_result <- function(object, ...) {
  head(object$result, ...)
}

#' head
#' @description Inspect crossdome result object
#' @noRd
View.xr_result <- function(object, ...) {
  View(object$result, ...)
}
