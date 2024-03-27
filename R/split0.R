#' Splits into a list
#'
#' Splits a vector (which is casted as a factor) by its elements returning a
#' named list with the indices for each unique element of the vector.
#'
#' @param x A vector to split.
#'
#' @return A named list according to the unique elements of `x` with the
#' integer indices of those given elements.
#'
#' @export
#' @author Andrew E Jaffe
#'
#' @seealso [splitit][rafalib::splitit]
#' @import rafalib
#'
#' @examples
#'
#' split0(letters[seq_len(3)])
#'
#' ## With some repeated info
#' set.seed(20161005)
#' abc <- sample(letters[seq_len(3)], 9, replace = TRUE)
#' split0(abc)
split0 <- function(x) {
    splitit(factor(x, levels = unique(x)))
}
