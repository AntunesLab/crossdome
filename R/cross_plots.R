#' cross_pairwise_plot
#'
#' @param query Description
#' @param subject Description
#'
#' @return Description
#'
#' @import pheatmap
#' @importFrom pheatmap pheatmap
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

  heatmap <- pheatmap::pheatmap(
    pairwise_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE
  )

  return(heatmap)
}
