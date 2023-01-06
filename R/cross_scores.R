#' Crossdome relatedness score
#'
#' @description Calculates relatedness score between peptides
#'
#' @param query A peptide target. Only 9-mers are supported.
#' @param subject Putative off-target candidate. Only 9-mers are supported.
#' @param position_weight A numeric vector derived from TCR hotspots
#'
#' @importFrom stats cor.test
#'
#' @return Returns the `list` holding Crossdome metrics
#' @export
#'
#' @examples
#' query <- 'EVDPIGHLY'
#' subject <- 'ESDPIVAQY'
#' cross_pair_summary(query = query, subject = subject)

cross_pair_summary <- function(query, subject, position_weight = rep(1, 9)) {

  query_vector <- .internal_checking_peptide(query)
  subject_vector <- .internal_checking_peptide(subject)

  if(length(query_vector) != length(subject_vector)) {
    quit("Please, the input sequence should have same length.")
  }

  n_mismatch <- length(query_vector) - base::sum(query_vector == subject_vector)
  n_positive <- base::sum(base::diag(BLOSUM80[query_vector, subject_vector]) > 0)

  query_components <- .internal_peptide_to_matrix(query_vector)
  subject_components <- .internal_peptide_to_matrix(subject_vector)

  if(length(position_weight) == length(query_vector)) {
    relatedness_score <- .internal_related_distance(query_components, subject_components, position_weight)
  } else {
    quit("Please, the position weight vector should match the peptides length.")
  }

  return(
    list(
      query = query,
      subject = subject,
      n_positive = n_positive,
      n_mismatch = n_mismatch,
      relatedness_score = relatedness_score
    )
  )

}
