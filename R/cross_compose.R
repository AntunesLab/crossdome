
#' Crossdome Compose
#'
#' @param query Description
#' @param subject Description
#' @param allele Description
#' @param position_weight Description
#'
#' @return Description
#' @export
#'
#' @examples
#' data('mage_off_targets')
#'
#' query <- 'EVDPIGHLY'
#' subject <- mage_off_targets$peptide_sequence
#' rank <- cross_compose(query = query, subject = subject)

cross_compose <- function(query, subject, allele, position_weight = NULL) {

  peptides <- subject$peptides

  if(length(peptides) > 1) {
    result <- lapply(peptides, function(off_target) {
      cross_pair_summary(query, off_target, position_weight)
    })
  } else {
    result <- cross_pair_summary(query, peptides, position_weight)
  }

  result <- do.call(rbind.data.frame, result)
  result <- result[order(result$diagonal_score, decreasing = T), ]

  result <- structure(
    list(
      position_weight = position_weight,
      allele = subject$allele,
      result = result,
      timestamp = utils::timestamp()
    ),
    class = "xr_result"
  )

  return(result)
}
