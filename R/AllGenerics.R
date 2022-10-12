#' cross_compose
#' @export
setGeneric('cross_compose', function(query, object, position_weight = rep(1, 9)) {
  standardGeneric("cross_compose")
})

#' cross_expression_matrix
#' @export
setGeneric('cross_expression_matrix', function(object, ...) {
  standardGeneric("cross_expression_matrix")
})

#' cross_expression_plot
#' @export
setGeneric('cross_expression_plot', function(object) {
  standardGeneric("cross_expression_plot")
})

#' cross_substitution_matrix
#' @export
setGeneric('cross_substitution_matrix', function(object, ...) {
  standardGeneric("cross_substitution_matrix")
})

#' cross_substitution_plot
#' @export
setGeneric('cross_substitution_plot', function(object) {
  standardGeneric("cross_substitution_plot")
})

#' cross_prediction_plot
#' @export
setGeneric('cross_prediction_plot', function(object) {
  standardGeneric("cross_prediction_plot")
})

#' cross_write
#' @export
setGeneric('cross_write', function(object, ...) {
  standardGeneric("cross_write")
})
