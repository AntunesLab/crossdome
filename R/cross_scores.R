#' Crossdome Summary
#'
#' @param query Description
#' @param subject Description
#' @param method Description
#'
#' @return
#' @export
#'
#' @examples
#' query <- 'EVDPIGHLY'
#' subject <- 'EVDPIGMLY'
#' cross_bp_summary(query = query, subject = subject)

cross_bp_summary <- function(query, subject, method = c("pearson", "kendall", "spearman")) {

  query_vector <- .internal_checking_epitope(query)
  subject_vector <- .internal_checking_epitope(subject)

  if(length(query_vector) != length(subject_vector)) {
    quit("Please, the input sequence should have same length.")
  }

  n_mismatch <- length(query_vector) - sum(query_vector == subject_vector)
  n_positive <- sum(base::diag(BLOSUM80[query_vector, subject_vector]) > 0)

  query_components <- .internal_epitope_to_matrix(query_vector)
  subject_components <- .internal_epitope_to_matrix(subject_vector)

  pairwise_matrix <- stats::cor(
    query_components, subject_components, method = method)

  pvalue <- cor.test(query_components, subject_components)$p.value
  diagonal_score <- sum(base::diag(pairwise_matrix)) / length(query)
  matrices_score <- .internal_matrix_metrics(query_components, subject_components)

  return(
    list(
      query = query,
      subject = subject,
      n_positive = n_positive,
      n_mismatch = n_mismatch,
      diagonal_score = diagonal_score,
      pvalue = pvalue,
      frobenius = matrices_score$frobenius,
      braun_score = matrices_score$braun_score
    )
  )
}
