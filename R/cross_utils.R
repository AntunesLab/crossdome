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
#' @param peptides Description
#' @param allele Description
#'
#' @return Description
#'
#' @import shiny
#' @importFrom shiny runApp
#' @export
#'
#' @examples
#' #' # Using default immunopeptidomics
#' background <- cross_universe(allele = 'HLA-A*01:01')
#'
#' # Using MAGE3A off-targets
#' data('mage_off_targets')
#'
#' mage_off_targets <- mage_off_targets$peptide_sequence
#' background <- cross_universe(off_targets = mage_off_targets, allele = "HLA-A*01:01")

cross_universe <- function(off_targets = NULL, allele) {

  hla_database <- crossdome::hla_database
  allele_list <- unique(hla_database$hla_allele)

  if(allele %in% allele_list) {
    background <- hla_database[
      hla_database$hla_allele == allele, 'peptide_sequence']
  }

  if(missing(off_targets)) {
    peptides <- background
  } else {
    peptides <- union(off_targets, background)
  }

  return(
    structure(
      list(
        allele = allele,
        peptides = peptides
        ),
      class = "xr_background"
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
#' @examples
#' peptides <- c("EVDPIGHLY", "ESDPIVAQY")
#' cross_target_expression(epitope = epitope)

cross_target_expression <- function(peptides) {

  if(missing(peptides)) {
    stop("Please, select a epitope sequence.")
  }

  peptide_annotation <- crossdome::peptide_annotation
  peptide_annotation <-
    peptide_annotation[match(peptides, peptide_annotation$peptide_sequence, nomatch = 0),]

  if(nrow(peptide_annotation) == 0) {
    stop("The provided epitopes are not included on the human proteome")
  } else {
    match_ratio <- setdiff(peptides, peptide_annotation$peptide_sequence)
    warning(
      paste0("Matching ration ", length(peptides) - length(match_ratio), " out of ", length(peptides), ". Not mapped epitopes. ",
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
