#' @name cross_pairwise_plot
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

  query_components <- .internal_peptide_to_matrix(query)
  subject_components <- .internal_peptide_to_matrix(subject)

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

#' @name cross_expression_heatmap
#'
#' @param object Description
#' @param pvalue_threshold Description
#'
#' @return Description
#'
#' @import ComplexHeatmap
#' @importFrom grid gpar unit
#'
#' @exportMethod cross_expression_heatmap
#'
#' @examples
#' if(FALSE) {
#'
#' result <- cross_expression_matrix(object = result)
#' cross_expression_heatmap(object = result)
#' }

setMethod(
  'cross_expression_heatmap', signature(object = "xrResult"),
  function(object) {

    if(nrow(object@result) == 0) {
      quit(
        paste0("The result object is empty. Check your threshold, pvalue_threshold ≤ ", pvalue_threshold)
      )
    }

    expression_matrix <-object@expression$data
    gene_donor_list <- expression_matrix$gene_donor

    expression_matrix <- as.matrix(
      expression_matrix[, 3:(ncol(expression_matrix)-1)]
    )

    rownames(expression_matrix) <- gene_donor_list
    expression_matrix <- expression_matrix[!duplicated(expression_matrix), ]

    p1 <- ComplexHeatmap::pheatmap(
      expression_matrix,
      scale = 'row',
      run_draw = FALSE
    )

    p2 <- ComplexHeatmap::rowAnnotation(
      Genes = ComplexHeatmap::anno_text(rownames(expression_matrix)),
      `mRNA expression` = ComplexHeatmap::anno_barplot(
        sqrt(rowSums(expression_matrix)),
        gp = grid::gpar(col = "black", fill = "#FFE200"),
        width = grid::unit(4, "cm"),
        border = FALSE,
        show_annotation_name = TRUE)
    )

    heatmap <- ComplexHeatmap::draw(p1 + p2)
    return(heatmap)
  }

)

#' @name cross_substitution_plot
#'
#' @param object Description
#' @param top Description
#'
#' @return Description
#'
#' @import ggplot2
#' @import patchwork
#' @importFrom universalmotif create_motif view_motifs
#' @importFrom Biostrings AAStringSet
#'
#' @exportMethod cross_substitution_plot
#'
#' @examples
#' if(FALSE) {
#'  cross_substitution_plot(object = result)
#' }


setMethod('cross_substitution_plot', signature(object = "xrResult"),
 function(object) {

   substitution <- object@analysis$substitution

   if(nrow(substitution) == 0) {
     quit(
       paste0("The substitution analysis not found. Consider to run `cross_substitution_matrix` method.")
     )
   }


   p1 <- universalmotif::view_motifs(prob_pos, use.type = 'PPM', sort.positions = T) +
     labs(y = "Probability") +
     theme(
       axis.text.x = element_text(size = 18, color = 'black'),
       axis.text.y = element_text(size = 14, color = 'black'),
       axis.title.y = element_text(size = 12)
     )

   substitution_pivoted <- .internal_prettify(substitution)
   substitution_pivoted$interval <- cut(
     substitution_pivoted$ppm,
     breaks = c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.85, 0.95, 1),
     include.lowest = TRUE
     )

   x_axis_sequence <- .internal_checking_peptide(object@query)
   x_axis_sequence <- setNames(1:9, x_axis_sequence)

   p2 <- ggplot(prob_pivoted, aes(x = position, y = reorder(aminoacid, -aa_idx), fill = ppm)) +
     geom_tile(colour = "grey", size = 0.45) +
     labs(x = NULL, y = 'Substitution', fill = "Probability") +
     scale_y_discrete(expand = c(0, 0)) +
     scale_x_discrete(expand = c(0,0), labels = x_axis_sequence, position = 'top') +
     theme_bw() +
     theme(
       axis.ticks = element_line(size = 0.4),
       axis.text = element_text(size = 12, color = 'black')
     )

   plot_seq <- p2 / p1 + plot_layout(heights = c(4, 1))
   return(plot_seq)

 }
)
