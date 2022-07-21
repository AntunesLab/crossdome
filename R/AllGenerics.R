#' cross_compose
#' @export
setGeneric('cross_compose', function(query, background, position_weight = rep(1, 9)) {
  standardGeneric("cross_compose")
})

#' cross_write
#' @export
setGeneric('cross_write', function(object, ...) {
  standardGeneric("cross_write")
})

#' #' cross_target_expression
#' #' @export
#' setGeneric('cross_target_expression', function(query, background, position_weight = rep(1, 9)) {
#'   standardGeneric("cross_compose")
#' })
