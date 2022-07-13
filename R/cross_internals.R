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

#' .internal_matrix_correlation
#' @description Returns Frobenius and Braun et al score
#' @noRd

.internal_matrix_correlation <- function(
  query_components,
  subject_components,
  method = c("pearson", "kendall", "spearman")
  ) {

  pairwise_matrix <- stats::cor(
    query_components, subject_components, method = method)
  pvalue <- cor.test(query_components, subject_components)$p.value
  diagonal_score <- sum(base::diag(pairwise_matrix)) / length(query_components)

  return(
    list(
      pairwise_matrix = pairwise_matrix,
      pvalue = pvalue,
      diagonal_score = diagonal_score)
    )
}

#' .internal_matrix_metrics
#' @description Returns Frobenius and Braun et al score
#' @noRd

.internal_matrix_metrics <- function(query_components, subject_components, position_weight = NULL) {

  product_components <- (query_components - subject_components)**2
  frobenius <- sum(product_components)
  braun_score <- sqrt(apply(product_components, 2, sum))

  if(!is.null(position_weight)) {
    position_weight_sqrt <- sqrt(position_weight)
    braun_score <- braun_score * position_weight_sqrt
  }

  braun_score <- sum(braun_score) / ncol(query_components)

  return(list(frobenius = frobenius, braun_score = braun_score))
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
