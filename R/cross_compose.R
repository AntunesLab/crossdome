
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
#' @import DT
#' @importFrom DT datatable
#'
#' @examples
#' data('mage_off_targets')
#'
#' query <- 'EVDPIGHLY'
#' subject <- mage_off_targets$peptide_sequence
#' rank <- cross_compose(query = query, subject = subject, widgets = FALSE)

cross_compose <- function(query, subject, allele, position_weight, widgets = FALSE) {

  if(length(subject) > 1) {
    result <- lapply(subject, function(off_target) {
      cross_bp_summary(query, off_target, position_weight)
    })
  } else {
    result <- cross_bp_summary(query, subject, position_weight)
  }

  result <- do.call(rbind.data.frame, result)
  result <- result[order(result$diagonal_score, decreasing = T), ]

  if(widgets) {
    cross_report_table(result)
  }

  return(result)
}
