#' @name cross_compose
#' @docType methods
#'
#' @param query Description
#' @param background Description
#' @param position_weight Description
#'
#' @return Description
#' @importFrom methods new
#' @exportMethod cross_compose
#'
#' @examples
#' data('mage_off_targets')
#'
#' mage_off_targets <- mage_off_targets$peptide_sequence
#' background <- cross_universe(off_targets = mage_off_targets, allele = "HLA-A*01:01")
#' result <- cross_compose(query = query, background = background)

setMethod(
  'cross_compose', signature(background = "xrBackground"),
  function(query, background, position_weight = rep(1, 9)) {

    peptides <- background@peptides

    if(length(peptides) > 1) {
      result <- lapply(peptides, function(off_target) {
        cross_pair_summary(query, off_target, position_weight)
      })
    } else {
      result <- cross_pair_summary(query, peptides, position_weight)
    }

    result <- do.call(rbind.data.frame, result)
    result <- result[order(result$relatedness_score, decreasing = F),]

    result <- new('xrResult',
                  result = result,
                  allele = background@allele,
                  position_weight = position_weight,
                  timestamp = utils::timestamp()
    )

    return(result)
  }
)
