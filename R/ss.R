#' Wrapper for string split and sapply
#'
#' Given a pattern, it splits and returns the slot of interest.
#'
#' @param x Passed to argument `x` of [strsplit][base::strsplit].
#' @param pattern Passed to argument `split` of [strsplit][base::strsplit].
#' @param slot An integer specifying which element of the resulting list to
#' return.
#' @param ... Additional arguments passed to [strsplit][base::strsplit].
#'
#' @return A vector with the information extracted from `x`.
#'
#' @export
#' @author Andrew E Jaffe
#'
#' @examples
#'
#' ## Some example info with two variables separated by a _
#' info <- c('DLPFC_polyA', 'DLPFC_RiboZero')
#' ss(info, '_', 1)
#' ss(info, '_', 2)
#'
#' ## Example where we use the ... argument
#' ss('a.b.c', '.', slot = 1, fixed = TRUE)
#'

ss <- function(x, pattern, slot = 1, ...) {
    sapply(strsplit(x = x, split = pattern, ...), '[', slot)
}
