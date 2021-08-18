#' Biochemistry Score
#'
#' @param query
#' @param subject
#'
#' @return
#' @export
#'
#' @examples
biochemistry_score <- function(query, subject) {
  query <- .internal_peptide_to_matrix(query)
  subject <- .internal_peptide_to_matrix(subject)
}
