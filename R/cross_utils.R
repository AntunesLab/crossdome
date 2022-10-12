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
  return(.internal_peptide_to_matrix(epitope))
}

#' cross_universe
#'
#' @param peptides Description
#' @param allele Description
#'
#' @return Description
#' @export
#' @examples
#' # Using default immunopeptidomics
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
#' @param object Description
#' @param pvalue_threshold Description
#'
#' @return Description
#'
#' @exportMethod cross_expression_matrix
#'
#' @examples
#'
#' result <- crossdome::mage_result
#' cross_expression_matrix(object = result)

setMethod(
  'cross_expression_matrix', signature(object = "xrResult"),
  function(object, pvalue_threshold = 0.01) {

    provisional <- object@result[
      object@result$pvalue <= pvalue_threshold, ]

    peptides <- provisional$subject
    peptide_annotation <- crossdome::peptide_annotation

    peptide_annotation <-
      peptide_annotation[match(peptides, peptide_annotation$peptide_sequence, nomatch = 0),]

    if(nrow(peptide_annotation) == 0) {
      stop("The provided epitopes are not included on the human proteome.")
    } else {
      match_ratio <- setdiff(peptides, peptide_annotation$peptide_sequence)
      warning(
        paste0("Matching ration ", length(peptides) - length(match_ratio), " out of ", length(peptides), ". Not mapped epitopes. ",
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
#' @exportMethod cross_substitution_matrix
#'
#' @examples
#' if(FALSE) {
#'  result <- cross_substitution_matrix(object = result)
#' }

setMethod('cross_substitution_matrix', signature(object = "xrResult"),
    function(object, pvalue_threshold = 0.01) {

      provisional <- object@result[
        object@result$pvalue <= pvalue_threshold, ]

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

#' show
#' @name show
#'
#' @param object Description
#' @docType method
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
#' @param object crossdome xrResult object
#' @param file 	path to output file.
#' @docType method
#' @exportMethod show
#' @importFrom utils write.table

setMethod('cross_write', signature(object = 'xrResult'),
          function(object, file = "", append = FALSE, quote = FALSE, sep = "\t",
                   eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE,
                   qmethod = c("escape", "double"), fileEncoding = "") {
            write.table(object@result, file = file, quote = quote, sep = sep, row.names = row.names)
          }
)


# cross_custom_prediction
# cross_custom_annotation
