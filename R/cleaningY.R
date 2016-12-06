#' Regress out SVs or PCs
#'
#' Regress out surrogate variables or principal components
#'
#' @param y The outcome
#' @param mod A model matrix
#' @param P The number of SVs or PCs to regress out.
#'
#' @return An object of the same type as `y` with the SVs/PCs regressed out.
#'
#' @export
#' @author Andrew E Jaffe
#'
#' @examples
#'
#' 

cleaningY <- function(y, mod, P = ncol(mod)) {
    stopifnot(P <= ncol(mod))
    Hat <- solve(t(mod) %*% mod) %*% t(mod)
    beta <- (Hat %*% t(y))
    cleany <- y - t(as.matrix(mod[, -c(seq_len(P))]) %*% beta[-seq_len(P), ])
    return(cleany)
}
