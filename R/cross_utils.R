#' cross_pairwise_plot
#'
#' @param query Description
#' @param subject Description
#'
#' @return
#'
#' @import ComplexHeatmap
#' @importFrom ComplexHeatmap pheatmap
#' @export
#'
#' @examples
#' query <- 'EVDPIGHLY'
#' subject <- 'EVDPIGMLY'
#' cross_pairwise_plot(query = query, subject = subject)

cross_pairwise_plot <- function(query, subject) {

  query_components <- .internal_epitope_to_matrix(query)
  subject_components <- .internal_epitope_to_matrix(subject)

  matrices_correlation <- .internal_matrix_correlation(
    query_components, subject_components)
  pairwise_matrix <- matrices_correlation$pairwise_matrix

    heatmap <- ComplexHeatmap::pheatmap(
      pairwise_matrix,
      cluster_rows = FALSE,
      cluster_cols = FALSE
    )

  return(heatmap)
}

#' cross_report_table
#'
#' @param result Description
#'
#' @return
#' @export
#'
#' @examples
#'

cross_report_table <- function(result) {
  return(datatable(result))
}

#' cross_epitope_properties
#'
#' @param epitope Description
#'
#' @return
#' @export
#'
#' @examples
#' cross_epitope_properties('EVDPIGHLY')

cross_epitope_properties <- function(epitope) {
  return(.internal_epitope_to_matrix(epitope))
}

#' cross_epitope_properties
#'
#' @param epitope Description
#'
#' @return
#'
#' @import shiny
#' @importFrom shiny runApp
#' @export
#'
#' @examples
#' cross_epitope_properties('EVDPIGHLY')

cross_browser <- function() {
  app_directory <- system.file("cross_browser", package = "crossdome")
  if (app_directory == "") {
    stop("Could not find example directory. Try re-installing `crossdome`.", call. = FALSE)
  }

  shiny::runApp(app_directory, display.mode = "normal")
}
