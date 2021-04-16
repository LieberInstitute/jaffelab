#' Calculate the odds ratio
#'
#' Calculate the odds ratio from a 2 by 2 table.
#'
#' @param x A two by two matrix.
#'
#' @return The odds ratio calculated from `x`.
#'
#' @export
#' @author Andrew E Jaffe
#'
#' @examples
#'
#' getOR(matrix(1:4, ncol = 2))
getOR <- function(x) {
    stopifnot(nrow(x) == 2)
    stopifnot(ncol(x) == 2)
    x[1, 1] / x[2, 1] / x[1, 2] * x[2, 2]
}
