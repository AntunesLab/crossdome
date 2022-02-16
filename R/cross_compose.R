
#' Crossdome Compose
#'
#' @param query Description
#' @param subject Description
#' @param allele Description
#' @param position_weight Description
#' @param widgets Description
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

  if(length(subject) > 1) {
    result <- lapply(subject, function(off_target) {
      cross_bp_summary(query, off_target, position_weight)
    })
  } else {
    result <- cross_bp_summary(query, subject, position_weight)
  }

  result <- do.call(rbind.data.frame, result)
  result <- result[order(result$diagonal_score, decreasing = T), ]

  return(result)
}
