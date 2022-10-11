#' .internal_checking_peptide
#' @importFrom Biostrings AA_STANDARD
#' @noRd

.internal_checking_peptide <- function(peptide) {

  peptide <- base::strsplit(peptide, split = "")[[1]]

  if(length(peptide) != 9) {
    quit("Please, Crossdome only supports 9-mer peptides.")
  }

  if(!all(peptide %in% Biostrings::AA_STANDARD)) {
    quit("Please, check your input sequence. Only standard aminoacids are supported.")
  }

  return(peptide)
}

#' .internal_peptide_to_matrix
#' @importFrom Biostrings AA_STANDARD
#' @importFrom stats cor
#' @noRd

.internal_peptide_to_matrix <- function(peptide) {

  if(length(peptide) <= 1) {
    peptide <- .internal_checking_peptide(peptide)
  }

  peptide_to_matrix <- MDS_COMPONENTS[, peptide]
  return(peptide_to_matrix)
}

#' .internal_related_distance
#' @noRd

.internal_related_distance <- function(query_components, subject_components, position_weight = NULL) {

  product_components <- (query_components - subject_components) ** 2
  relatedness_score <- sqrt(apply(product_components, 2, sum))

  if(!is.null(position_weight)) {
    relatedness_score <- relatedness_score * sqrt(position_weight)
  }

  relatedness_score <- sum(relatedness_score) / ncol(query_components)
  return(relatedness_score)
}

#' .internal_percentil_rank
#' @noRd

.internal_percentile_rank <- function(relatedness) {

  rank <- ((1:length(relatedness) - 1) / (sum(!is.na(relatedness)) -  1)) * 100

  return(rank)
}


#' .internal_prettify
#' @noRd

.internal_prettify <- function(ppm_matrix) {

  ppm_data <- sapply(1:length(rownames(ppm_matrix)), function(x) {
    data.frame(
      aa_idx = x,
      position = 1:9,
      aminoacid = rownames(ppm_matrix)[x],
      ppm = ppm_matrix[x, ]
    )
  }, simplify = FALSE)

  ppm_data <- do.call('rbind', ppm_data)

  return(ppm_data)

}
