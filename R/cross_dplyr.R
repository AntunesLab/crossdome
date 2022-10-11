#' @method select xrResult
#' @importFrom dplyr select
#' @importFrom dplyr quos
#' @export

select.xrResult <- function(.data, ..., .preserve = FALSE) {

  dots <- dplyr::quos(...)
  .data@result <- .data@result %>%
    dplyr::select(!!!dots)

  return(.data)

}

#' @method filter xrResult
#' @importFrom dplyr filter
#' @importFrom dplyr quos
#' @export

filter.xrResult <- function(.data, ..., .preserve = FALSE) {

  dots <- dplyr::quos(...)
  .data@result <- .data@result %>%
    dplyr::filter(!!!dots, .preserve = .preserve)

  return(.data)

}

#' @method mutate xrResult
#' @importFrom dplyr mutate
#' @importFrom dplyr quos
#' @export

mutate.xrResult <- function(.data, ...) {

  dots <- dplyr::quos(...)
  .data@result <- .data@result %>%
    dplyr::mutate(!!!dots)

  return(.data)

}
