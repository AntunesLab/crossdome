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
  'cross_compose', signature(object = "xrBackground"),
  function(query, object, position_weight = rep(1, 9)) {

    peptides <- object@peptides

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
    result$hla_allele <- object@allele

    result <- result[order(result$relatedness_score, decreasing = F), ]
    result$percentile_rank <- .internal_percentile_rank(result$relatedness_score)
    result$index <- 1:nrow(result)

    result <- result[,
          c('index', 'query', 'subject', 'n_positive', 'n_mismatch', 'relatedness_score',
            'zscore', 'pvalue', 'hla_allele', 'percentile_rank')]

    result <- new('xrResult',
                  query = query,
                  result = result,
                  allele = object@allele,
                  position_weight = position_weight,
                  timestamp = utils::timestamp()
    )

    return(result)

  }
)
