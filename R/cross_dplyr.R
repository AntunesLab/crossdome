#' @method select xrResult
#' @importFrom dplyr select
#' @improtFrom dplyr quos
#' @export

select.xrResult <- function(.data, ..., .preserve = FALSE) {

  dots <- dplyr::quos(...)
  .data@result <- .data@result %>%
    dplyr::select(!!!dots)

  return(.data)

}

#' @method filter xrResult
#' @importFrom dplyr filter
#' @export

filter.xrResult <- function(.data, ..., .preserve = FALSE) {

  dots <- dplyr::quos(...)
  .data@result <- .data@result %>%
    dplyr::filter(!!!dots, .preserve = .preserve)

  return(.data)

}

#' @method mutate xrResult
#' @importFrom dplyr mutate
#' @export

mutate.xrResult <- function(.data, ...) {

  dots <- dplyr::quos(...)
  .data@result <- .data@result %>%
    dplyr::mutate(!!!dots)

  return(.data)

}
