
#' Title
#'
#' @param query Description
#' @param subject Description
#' @param hla_allele Description
#' @param tcr_weights Description
#'
#' @return
#' @export
#'
#' @examples
#'

cross_compose <- function(query,
                          universe = HLA_LIST,
                          tcr_weights,
                          widget = TRUE) {


  results <- sapply(universe, cross_bp_summary, )

  if(widgets) {
    cross_report_table(results)
  }
  return(results)

}
