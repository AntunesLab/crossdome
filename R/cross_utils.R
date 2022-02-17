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

#' cross_epitope_properties
#'
#' @param epitope Description
#'
#' @return Description
#' @export
#'
#' @examples
#' cross_epitope_properties('EVDPIGHLY')

cross_epitope_properties <- function(epitope) {
  return(.internal_epitope_to_matrix(epitope))
}

#' cross_universe
#'
#' @param subject Description
#' @param allele Description
#'
#' @return Description
#'
#' @import shiny
#' @importFrom shiny runApp
#' @export
#'
#' @examples
#' subject <- mage_off_targets$peptide_sequence
#' cross_universe(subject, allele = 'HLA*A-01:01')

cross_universe <- function(subject, allele) {

  hla_database <- crossdome::hla_database
  allele_list <- unique(hla_database$hla_allele)

  if(allele %in% allele_list) {
    background <- hla_database[
      hla_database$hla_allele == allele, 'peptide_sequence']
  }

  if(missing(subject)) {
    subject <- background
  } else {
    subject <- union(subject, background)
  }

  return(
    list(
      allele = allele,
      subject = subject
    )
  )
}

#' cross_expression
#'
#' @param epitope Description
#' @param
#'
#' @return Description
#'
#' @export

cross_expression <- function(epitope) {
  peptide_annotation <- crossdome::peptide_annotation
}


#' cross_browser
#'
#' @return Description
#'
#' @import shiny
#' @importFrom DT dataTableOutput
#' @importFrom DT renderDataTable
#' @importFrom shiny runApp
#' @export

cross_browser <- function() {
  app_directory <- system.file("cross_browser", package = "crossdome")
  if (app_directory == "") {
    stop("Could not find example directory. Try re-installing `crossdome`.", call. = FALSE)
  }
  shiny::runApp(app_directory, display.mode = "normal")
}
