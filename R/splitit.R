#' Splits into a list
#'
#' Splits a vector by its elements returning a named list with the indices
#' for each unique element of the vector.
#'
#' @param x A vector to split.
#'
#' @return A named list according to the unique elements of \code{x} with the
#' integer indices of those given elements.
#'
#' @export
#' @author Andrew E Jaffe
#'
#' @seealso \link{split0}
#'
#' @examples
#'
#' splitit(letters[1:3])
#'
#' ## With some repeated info
#' set.seed(20161005)
#' abc <- sample(letters[1:3], 9, replace = TRUE)
#' splitit(abc)
#'

splitit <- function(x) {
    split(seq(along = x), x) 
}
