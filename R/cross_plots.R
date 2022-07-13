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
#' subject <- 'ESDPIVAQY'
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

#' cross_expression_heatmap
#'
#' @param result Description
#' @param epitope Description
#' @param top Description
#'
#' @return Description
#'
#' @importFrom ComplexHeatmap pheatmap
#' @importFrom ComplexHeatmap rowAnnotation
#' @importFrom ComplexHeatmap draw
#' @importFrom ComplexHeatmap anno_text
#' @export
#'
#' @examples
#' epitope <- c("EVDPIGHLY", "ESDPIVAQY")
#' cross_expression_heatmap(epitope = epitope)

cross_expression_heatmap <- function(result, epitope, top = 50) {

  if(!missing(result)) {
    epitope <- head(result[, 'subject'], n = top)
  }

  if(!missing(epitope)) {
    expression_df <- cross_target_expression(epitope)
  }

  expression_matrix <- as.matrix(
    expression_df[, 3:(ncol(expression_df)-1)]
  )
  rownames(expression_matrix) <- expression_df$gene_donor

  p1 <- ComplexHeatmap::pheatmap(
    expression_matrix,
    scale = 'row',
    run_draw = FALSE
  )

  p2 <- ComplexHeatmap::rowAnnotation(
    Genes = anno_text(rownames(expression_matrix)),
    `mRNA expression` = anno_barplot(
      sqrt(rowSums(expression_matrix)),
      gp = gpar(col = "black", fill = "#FFE200"),
      width = unit(4, "cm"),
      border = FALSE,
      show_annotation_name = TRUE)
  )

  heatmap <- ComplexHeatmap::draw(p1 + p2)
  return(heatmap)
}
