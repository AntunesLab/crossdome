
#' Crossdome Compose
#'
#' @param query Description
#' @param subject Description
#' @param widgets Description
#'
#' @return
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

cross_compose <- function(query, subject, widgets = FALSE) {

  if(length(subject) > 1) {
    result <- lapply(subject, function(off_target) {
      cross_bp_summary(query, off_target)
    })
  } else {
    cross_bp_summary(query, subject)
  }

  result <- do.call(rbind.data.frame, result)
  result <- result[order(result$diagonal_score, decreasing = T), ]

  if(widgets) {
    cross_report_table(result)
  }

  return(result)
}
