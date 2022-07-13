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
#' cross_universe(subject, allele = 'HLA-A*01:01')

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
    structure(
      list(allele = allele,
           peptides = subject),
      class = "xr_universe"
    )
  )
}

#' cross_target_expression
#'
#' @param epitope Description
#'
#' @return Description
#'
#' @export
#'
#' @example
#' epitope <- c("EVDPIGHLY", "ESDPIVAQY")
#' cross_target_expression(epitope = epitope)

cross_target_expression <- function(epitope) {

  if(missing(epitope)) {
    stop("Please, select a epitope sequence.")
  }

  peptide_annotation <- crossdome::peptide_annotation
  peptide_annotation <-
    peptide_annotation[match(epitope, peptide_annotation$peptide_sequence, nomatch = 0),]

  if(nrow(peptide_annotation) == 0) {
    stop("The provided epitopes are not included on the human proteome")
  } else {
    match_ratio <- setdiff(epitope, peptide_annotation$peptide_sequence)
    warning(
      paste0("Matching ration ", length(epitope) - length(match_ratio), " out of ", length(epitope), ". Not mapped epitopes. ",
             paste0(match_ratio, collapse = ",")
      )
    )
  }

  target_expression <- crossdome::gtex_database
  target_expression <- merge(
    target_expression,
    peptide_annotation,
    by = c('ensembl_id', 'gene_donor')
  )

  return(target_expression)
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
