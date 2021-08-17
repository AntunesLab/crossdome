#' Title
#' @description Returns MDS components for each amine acid
#' @noRd
.internal_checking_peptide <- function(peptide) {
  peptide <- base::strsplit(peptide, split = "")[[1]]
  if(!all(peptide %in% Biostrings::AA_STANDARD)) {
    quit("Please, check your input sequence.")
  }
  return(peptide)
}

#' Title
#' @description Returns MDS components for each amine acid
#' @importFrom Biostrings AA_STANDARD
#' @noRd

.internal_peptide_to_matrix <- function(peptide) {
  peptide <- .internal_checking_peptide(peptide)
  peptide_to_matrix <- MDS_COMPONENTS[, peptide]
  return(peptide_to_matrix)
}

#' Biochemistry Correlation
#'
#' @param query Query sequence
#' @param subject Subject sequence
#' @param method Spearman correlation as default
#'
#' @return
#' @export
#'
#' @examples
#' query <- 'EVDPIGHLY'
#' subject <- 'EVDPIGMLY'
#' biochemistry_correlation(query = query, subject = subject)

biochemistry_correlation <- function(query, subject, method = "spearman") {
  query <- .internal_peptide_to_matrix(query)
  subject <- .internal_peptide_to_matrix(subject)

  if(length(query) != length(subject)) {
    quit("Please, check your input sequence.")
  }

  pairwise_matrix <- stats::cor(query, subject, method = method)
  diagonal_score <- sum(base::diag(pairwise_matrix)) / length(query)

  cross_object <- structure(
    list(
      pairwise_matrix = pairwise_matrix,
      diagonal_score = diagonal_score
    ),
    class = 'cross_object'
  )

  return(cross_object)
}

#' Peptide Alignment Stats
#'
#' @param query
#' @param subject
#'
#' @return
#' @importFrom Biostrings pairwiseAlignment
#' @export
#'
#' @examples
#' query <- 'EVDPIGHLY'
#' subject <- 'EVDPIGMLY'
#' peptide_alignment_stats(query = query, subject = subject)

peptide_alignment_stats <- function(query, subject) {
  query <- .internal_checking_peptide(query)
  subject <- .internal_checking_peptide(subject)

  if(length(query) != length(subject)) {
    quit("Please, check your input sequence.")
  }

  n_mismatch <- length(query) - sum(query == subject)
  n_positive <- sum(base::diag(BLOSUM80[query, subject]) > 0)

  cross_alg_stats <- structure(
    list(
      n_mismatch = n_mismatch,
      n_positive = n_positive),
    class = "cross_alg_stats"
  )

  return(cross_alg_stats)
}
