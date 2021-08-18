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
  peptide <- .internal_checking_peptide(peptide)
  peptide_to_matrix <- MDS_COMPONENTS[, peptide]
  return(peptide_to_matrix)
}
