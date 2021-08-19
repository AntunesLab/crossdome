#' .internal_checking_peptide
#' @description Returns MDS components for each amine acid
#' @importFrom Biostrings AA_STANDARD
#' @noRd
.internal_checking_peptide <- function(peptide) {
  peptide <- base::strsplit(peptide, split = "")[[1]]
  if(!all(peptide %in% Biostrings::AA_STANDARD)) {
    quit("Please, check your input sequence.")
  }
  return(peptide)
}

#' .internal_peptide_to_matrix
#' @description Returns MDS components for each amine acid
#' @importFrom Biostrings AA_STANDARD
#' @noRd

.internal_peptide_to_matrix <- function(peptide) {
  peptide_to_matrix <- MDS_COMPONENTS[, peptide]
  return(peptide_to_matrix)
}


#' .internal_matrix_metrics
#' @description Returns Frobenius and Braun et al score
#' @noRd

.internal_matrix_metrics <- function(query_components, subject_components) {
  product_components <- (query_components - subject_components)**2

  frobenius <- sum(product_components)
  braun_score <- sum(
    sqrt(apply(product_components, 2, sum))
  )

  return(list(frobenius = frobenius, braun_score = braun_score))

}
