#' @name cross_compose
#' @title Performing Crossdome screening
#' @export
setGeneric('cross_compose', function(query, background, position_weight = rep(1, 9)) {
  standardGeneric("cross_compose")
})

#' @name cross_expression_matrix
#' @title Collecting gene-donor expression data
#' @export
#' @usage NULL
setGeneric('cross_expression_matrix', function(object, ...) {
  standardGeneric("cross_expression_matrix")
})

#' @name cross_expression_plot
#' @title Expression heatmap
#' @export
setGeneric('cross_expression_plot', function(object) {
  standardGeneric("cross_expression_plot")
})

#' @name cross_substitution_matrix
#' @title Substitution matrix analysis
#' @export
setGeneric('cross_substitution_matrix', function(object, ...) {
  standardGeneric("cross_substitution_matrix")
})

#' @name cross_substitution_plot
#' @title Substitution heatmap
#' @export
setGeneric('cross_substitution_plot', function(object) {
  standardGeneric("cross_substitution_plot")
})

#' @name cross_tissues_plot
#' @title Tissue-specificity barplot
#' @export
setGeneric('cross_tissues_plot', function(object) {
  standardGeneric("cross_tissues_plot")
})

#' @name cross_prediction_plot
#' @title Prediction dot plots
#' @export
setGeneric('cross_prediction_plot', function(object) {
  standardGeneric("cross_prediction_plot")
})

#' @name cross_write
#' @title Data Output
#' @export
setGeneric('cross_write', function(object, ...) {
  standardGeneric("cross_write")
})
