#' @name cross_pairwise_plot
#' @title Correlation plot
#'
#' @description Plot a Correlation plot based on two biochemical profiles
#'
#' @param query Peptide target. Only 9-mers are supported.
#' @param subject Putative off-target candidate. Only 9-mers are supported.
#'
#' @return Returns a \strong{corrplot} object
#'
#' @import corrplot
#' @importFrom corrplot corrplot
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
    query_components, subject_components
  )

  corplot <- corrplot::corrplot(
    matrices_correlation,
    title = paste0("\n\n", query, "-", subject),
    method = "color",
    tl.col = "black",
    shade.col	= 2,
    tl.srt = 45,
    type = "lower"
  )

  return(corplot)
}

#' @name cross_expression_plot
#'
#' @description Plot a heatmap presenting the gene donor expression profile
#'
#' @param object Depends on xrResult object. Run \code{\link{cross_compose}} function.
#'
#' @return Returns a combined \strong{ComplexHeatmap} object
#'
#' @importFrom ComplexHeatmap pheatmap rowAnnotation anno_text anno_barplot draw
#' @importFrom grid gpar unit
#'
#' @exportMethod cross_expression_plot
#'
#' @examples
#' \dontrun{
#'
#' result <- cross_expression_matrix(object = result)
#' cross_expression_plot(object = result)
#' }

setMethod('cross_expression_plot', signature(object = "xrResult"),
          function(object) {

            if(nrow(object@result) == 0) {
              quit(
                paste0("The result object is empty. Check your threshold, pvalue_threshold ≤ ", pvalue_threshold)
              )
            }

            expression_matrix <-object@expression$data
            gene_donor_list <- expression_matrix$gene_donor

            expression_matrix <- as.matrix(
              expression_matrix[, 4:(ncol(expression_matrix)-3)]
            )

            rownames(expression_matrix) <- gene_donor_list
            expression_matrix <- expression_matrix[!duplicated(expression_matrix), ]

            p1 <- ComplexHeatmap::pheatmap(
              expression_matrix,
              name = "mRNA Z-score",
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
#' @description Plot a heatmap combined with seqlogo displaying amino acid substitutions
#'
#' @param object Depends on xrResult object. Run \code{\link{cross_compose}} function.
#'
#' @return Returns a heatmap-like figure showing substitution rate
#'
#' @import ggplot2
#' @import patchwork
#' @importFrom universalmotif create_motif view_motifs
#' @importFrom Biostrings AAStringSet
#' @importFrom stats reorder setNames
#'
#' @exportMethod cross_substitution_plot
#'
#' @examples
#' \dontrun{
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

            p1 <- universalmotif::view_motifs(substitution, use.type = 'PPM', sort.positions = T) +
              labs(y = "Probability") +
              theme(
                axis.text.x = element_text(size = 18, color = 'black'),
                axis.text.y = element_text(size = 14, color = 'black'),
                axis.title.y = element_text(size = 12)
              )

            substitution_pivoted <- .internal_prettify(substitution)

            x_axis_sequence <- .internal_checking_peptide(object@query)
            x_axis_sequence <- setNames(
              x_axis_sequence, 1:length(x_axis_sequence))

            query_sequence <- merge(
              substitution_pivoted,
              data.frame(
                position = names(x_axis_sequence),
                aminoacid = x_axis_sequence
              ),
              by = c('position', 'aminoacid')
            )

            p2 <- ggplot(substitution_pivoted, aes(x = as.factor(position), y = reorder(aminoacid, aa_idx), fill = ppm)) +
              geom_tile(colour = "grey", size = 0.45) +
              geom_rect(data = query_sequence, size = 1, fill = NA, colour = "black",
                        aes(xmin = position - 0.5, xmax = position + 0.5, ymin = aa_idx - 0.5, ymax = aa_idx + 0.5)) +
              labs(x = NULL, y = 'Substitution', fill = "Probability") +
              scale_y_discrete(expand = c(0, 0)) +
              scale_x_discrete(expand = c(0, 0), labels = x_axis_sequence, position = 'top') +
              scale_colour_manual(values = c('Yes' = 'black', 'No' = 'grey')) +
              scale_fill_distiller(
                type = 'div',
                palette = 'BuPu',
                direction = 1,
                limits = c(0, 1),
                breaks = seq(0, 1, by = 0.25)) +
              theme_bw() +
              theme(
                axis.ticks = element_line(size = 0.4),
                axis.text = element_text(size = 14, color = 'black'),
                axis.title.x = element_text(size = 14)
              )


            plot_seq <- p2 / p1 + plot_layout(heights = c(4, 1))
            return(plot_seq)
          }
)

#' @name cross_enrichment_plot
#'
#' @description Plot a bar plot summarizing the tissue-specificy groups
#'
#' @param object Depends on xrResult object. Run \code{\link{cross_compose}} function.
#'
#' @return Returns a ggplot object
#'
#' @import ggplot2
#' @importFrom stats aggregate reorder
#' @importFrom methods slotNames
#'
#' @exportMethod cross_prediction_plot
#'
#' @examples
#' \dontrun{
#'  cross_enrichment_plot(object = result)
#' }

setMethod('cross_enrichment_plot', signature(object = "xrResult"),
          function(object) {
            if(!'expression' %in% slotNames(object)) {
              quit(
                paste0("Expression slot not found! Please, run the cross_expression_matrix function.")
              )
            }

            specificity_group <- object@expression$data[, c('Group', 'tissues')]
            specificity_group <- aggregate(tissues ~ ., data = specificity_group, length)

            ggplot(specificity_group, aes(x = reorder(Group, -tissues), y = tissues)) +
              geom_col() +
              labs(x = NULL, y = '# of candidates') +
              theme_classic() +
              theme(
                axis.text = element_text(size = 12, colour = 'black')
              )
          }
)

#' @name cross_prediction_plot
#'
#' @description Plot a dot plot showing the immunogenic predictions
#'
#' @param object Depends on xrResult object. Run \code{\link{cross_compose}} function.
#'
#' @return Returns a ggplot object
#'
#' @import ggplot2
#' @import patchwork
#' @importFrom universalmotif create_motif view_motifs
#' @importFrom Biostrings AAStringSet
#'
#' @exportMethod cross_prediction_plot
#'
#' @examples
#' \dontrun{
#'  cross_prediction_plot(object = result)
#' }


setMethod('cross_prediction_plot', signature(object = "xrResult"),
          function(object) {
            print("Under construction...")
          }
)
