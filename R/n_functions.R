#' Get more than just the single max or single min value with n_max() and n_min()
#'
#' Returns the n maximum or n minimum values of a vector.
#'
#' @param vec A vector from which to select the max or min values.
#' @param n The number of values to return.
#'
#' @return The n max or min values of your vector.
#' @export
#' @author Emily E. Burke
#' @rdname n_functions
#'
#' @examples
#'
#' x <- rnorm(100)
#' n_max(x)
#' n_min(x)
n_max <- function(vec, n = 5) {
    return(head(sort(vec, decreasing = TRUE), n))
}


#' @rdname n_functions
#' @export
n_min <- function(vec, n = 5) {
    return(head(sort(vec, decreasing = FALSE), n))
}
