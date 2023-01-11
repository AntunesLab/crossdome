#' @name cross_peptide_properties
#' @title Converts a peptide to biochemical profile
#'
#' @param query A peptide target. Only 9-mers are supported.
#'
#' @return Returns a matrix related to the biochemical profile.
#' @export
#'
#' @examples
#' cross_peptide_properties('EVDPIGHLY')

cross_peptide_properties <- function(query) {
  return(.internal_peptide_to_matrix(query))
}

#' @name cross_background
#' @title Creating immunopeptidomics database
#'
#' @description Peptide database spanning eluted candidates (experimentally validated) and custom (user-defined).
#'
#' @param off_targets A list of off-target candidates. Only 9-mers are supported.
#' @param allele Input an MHC Class I allele. Please, check \code{\link{hla_database}}.
#'
#' @return Returns a \code{\link{xrBackground}} object
#' @export
#' @examples
#' \dontrun{
#'
#' # Listing MHC class I alleles
#' data('hla_database')
#' View(hla_database)
#'
#' # Using default immunopeptidomics
#' background <- cross_background(allele = 'HLA-A*01:01')
#'
#' # Using MAGE3A off-targets
#' data('mage_off_targets')
#'
#' mage_off_targets <- mage_off_targets$peptide_sequence
#' background <- cross_background(off_targets = mage_off_targets, allele = "HLA-A*01:01")
#'
#'}

cross_background <- function(off_targets = NULL, allele) {

  hla_database <- crossdome::hla_database
  allele_list <- unique(hla_database$hla_allele)

  if(allele %in% allele_list) {

    background <- hla_database[
      hla_database$hla_allele == allele, 'peptide_sequence']

  } else {

    quit(paste0("Crossdome does not supports the ", allele, ". For more details check the documentation."))

  }

  if(!missing(off_targets)) {

    background <- union(off_targets, background)

  }

  object <- new('xrBackground',
                  allele = allele,
                  peptides = background,
                  stats = list(
                    'off-target' = ifelse(is.null(off_targets), 0, length(off_targets)),
                    'database' = length(background)
                  )
  )

  return(object)

}

#' @name cross_expression_matrix
#' @docType methods
#'
#' @description Extracts gene donor mRNA expression based on CR candidates.
#'
#' @param object Depends on xrResult object. Run \code{\link{cross_compose}} function.
#' @param rank A numeric value to filter number of candidates
#' @param pvalue_threshold P-value threshold
#'
#' @return Return a matrix containing mRNA expression across healthy tissues.
#'
#' @importFrom utils head
#'
#' @exportMethod cross_expression_matrix
#'
#' @examples
#' \dontrun{
#' result <- crossdome::mage_result
#' cross_expression_matrix(object = result)
#' }

setMethod(
  'cross_expression_matrix', signature(object = "xrResult"),
  function(object, rank = NULL, pvalue_threshold = 0.01) {

    provisional <- object@result[
      object@result$pvalue <= pvalue_threshold, ]

    if(!is.null(rank)) {
      provisional <- head(provisional, n = rank)
    }

    peptides <- provisional$subject
    peptide_annotation <- crossdome::peptide_annotation

    peptide_annotation <-
      peptide_annotation[match(peptides, peptide_annotation$peptide_sequence, nomatch = 0),]

    if(nrow(peptide_annotation) == 0) {
      stop("The provided epitopes are not included on the human proteome.")
    } else {
      match_ratio <- setdiff(peptides, peptide_annotation$peptide_sequence)
      warning(
        paste0("Found ", length(peptides) - length(match_ratio), " peptides out of ", length(peptides), ". Unmapped peptides: ",
               paste0(match_ratio, collapse = ",")
        )
      )
    }

    target_expression <- crossdome::hpa_database
    target_expression <- merge(
      target_expression,
      peptide_annotation,
      by = c('ensembl_id', 'gene_donor')
    )

    target_expression <- target_expression[, c(1, 2, 43, 3:42)]

    object@expression = list(
      'data' = target_expression,
      'umapped' = match_ratio,
      'pvalue' = pvalue_threshold
    )

    return(object)

  }
)

#' @name cross_substitution_matrix
#'
#' @description Calculates position-specific substitution across cross-reactive candidates.
#'
#' @param object Depends on xrResult object. Run \code{\link{cross_compose}} function.
#' @param rank A numeric value to filter number of candidates
#' @param pvalue_threshold P-value threshold
#'
#' @return Returns a \code{matrix} with amino acid substitution probabilies.
#'
#' @importFrom universalmotif create_motif view_motifs
#' @importFrom Biostrings AAStringSet
#' @importFrom utils head
#'
#' @exportMethod cross_substitution_matrix
#'
#' @examples
#' \dontrun{
#'  result <- cross_substitution_matrix(object = result)
#' }

setMethod('cross_substitution_matrix', signature(object = "xrResult"),
    function(object, rank = NULL, pvalue_threshold = 0.01) {

      provisional <- object@result[
        object@result$pvalue <= pvalue_threshold, ]

      if(!is.null(rank)) {
        provisional <- head(provisional, n = rank)
      }

      if(nrow(provisional) == 0) {
        quit(
          paste0("The result object is empty. Check your threshold, pvalue_threshold ≤", pvalue_threshold)
        )
      }

      peptides <- Biostrings::AAStringSet(provisional$subject)
      prob_pos <- universalmotif::create_motif(peptides, type = 'PPM')

      object@analysis <- list(
        'substitution' = prob_pos@motif,
        'pvalue' = pvalue_threshold
      )

      return(object)

    }
)

#' @name cross_browser
#' @title Launching Crossdome web application
#'
#' @description Opens an interactive shiny application.
#'
#'
#' @importFrom DT dataTableOutput
#' @importFrom DT renderDataTable
#' @importFrom shiny runApp
#'
#' @export
#' @examples
#' \dontrun{
#'
#'  # Opening the shiny application
#'  cross_browser()
#'
#' }


cross_browser <- function() {
  app_directory <- system.file("cross_browser", package = "crossdome")
  if (app_directory == "") {
    stop("Could not find example directory. Try re-installing `crossdome`.", call. = FALSE)
  }
  shiny::runApp(app_directory, display.mode = "normal")
}

#' @name show
#' @title show
#'
#' @param object Depends on xrResult object. Run \code{\link{cross_compose}} function.
#'
#' @exportMethod show
#' @importFrom utils View

setMethod('show', signature(object = 'xrResult'),
          function(object) {
            return(
              utils::View(object@result, title = 'Result')
            )
          }
)

#' @name cross_write
#'
#' @description Export Crossdome object to a tsv file.
#'
#' @param object Depends on xrResult object. Run \code{\link{cross_compose}} function.
#' @param file File or connection to write to.
#'
#' @exportMethod show
#' @importFrom utils write.table

setMethod('cross_write', signature(object = 'xrResult'),
          function(object, file = "", append = FALSE, quote = FALSE, sep = "\t",
                   eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE,
                   qmethod = c("escape", "double"), fileEncoding = "") {
            write.table(object@result, file = file, quote = quote, sep = sep, row.names = row.names)
          }
)
