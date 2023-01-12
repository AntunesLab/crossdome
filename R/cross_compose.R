#' @name cross_compose
#' @docType methods
#'
#' @description Predicts relatedness among peptides in a given database. Low values are associated with cross-reactive candidates.
#'
#' @param query A peptide target. Only 9-mers are supported.
#' @param background Background object create with \code{\link{cross_background}}
#' @param position_weight A numeric vector derived from TCR hotspots
#'
#' @return Returns a \code{\link{xrResult}} object containing relatedness ranking
#' @importFrom methods new
#' @exportMethod cross_compose
#'
#' @examples
#' \dontrun{
#'
#' query <- 'EVDPIGHLY'
#' data('mage_off_targets')
#'
#' mage_off_targets <- mage_off_targets$peptide_sequence
#' background <- cross_background(off_targets = mage_off_targets, allele = "HLA-A*01:01")
#' result <- cross_compose(query = query, background = background)
#'
#' }

setMethod(
  'cross_compose', signature(background = "xrBackground"),
  function(query, background, position_weight = rep(1, 9)) {

    peptides <- background@peptides

    if(length(peptides) > 1) {
      result <- lapply(peptides, function(off_target) {
        crossdome::cross_pair_summary(query, off_target, position_weight)
      })
    } else {
      result <- crossdome::cross_pair_summary(query, peptides, position_weight)
    }

    result <- do.call('rbind.data.frame', result)

    result$zscore <- base::scale(result$relatedness_score, center = 31.97307, scale = 5.834203)[,1]
    result$pvalue <- stats::pnorm(result$zscore, lower.tail = TRUE)
    result$hla_allele <- background@allele

    result <- result[order(result$relatedness_score, decreasing = F), ]
    result$percentile_rank <- .internal_percentile_rank(result$relatedness_score)
    result$rank <- 1:nrow(result)

    result <- result[,
          c('rank', 'query', 'subject', 'n_positive', 'n_mismatch', 'relatedness_score',
            'zscore', 'pvalue', 'hla_allele', 'percentile_rank')]

    result <- new('xrResult',
                  query = query,
                  result = result,
                  allele = background@allele,
                  position_weight = position_weight,
                  timestamp = utils::timestamp()
    )

    return(result)

  }
)
