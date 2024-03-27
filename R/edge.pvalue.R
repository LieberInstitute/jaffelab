#' Function to compute empirical p-values
#'
#' Given a pattern, it splits and returns the slot of interest.
#'
#' @param stat vector of observed statistics.
#' @param stat0 matrix of null statistics, where rows are features and columns
#'  are iterations.
#' @param pool features are interchangeable, e.g. unlist all null statistics.
#'
#' @return A vector of empirical p-values.
#'
#' @export
#' @author Jeffrey T Leek, Leonardo Collado-Torres (examples)
#'
#' @examples
#'
#' ## Vector of statistics
#' z <- rnorm(100, mean = 2)
#'
#' ## Matrix with one row per element in the vector of statistics
#' ## and one column per permutation (here 100)
#' z0 <- matrix(rnorm(1e5), nrow = 100)
#'
#' ## Pool all the permutations (default)
#' edge.pvalue(z, z0, pool = TRUE)
#'
#' ## Or don't pull them
#' edge.pvalue(z, z0, pool = FALSE)
edge.pvalue <- function(stat, stat0, pool = TRUE) {
    err.func <- "edge.pvalue"
    m <- length(stat)
    if (pool == TRUE) {
        if (is.matrix(stat0)) {
            stat0 <- as.vector(stat0)
        }
        m0 <- length(stat0)
        v <- c(rep(TRUE, m), rep(FALSE, m0))
        v <- v[order(c(stat, stat0), decreasing = TRUE)]
        u <- seq_len(length(v))
        w <- seq_len(m)
        p <- (u[v == TRUE] - w) / m0
        p <- p[rank(-stat)]
        p <- pmax(p, 1 / m0)
    } else {
        if (is.vector(stat0)) {
            stop("stat0 must be a matrix.", call. = FALSE)
        }
        if (ncol(stat0) == m) {
            stat0 <- t(stat0)
        }
        if (nrow(stat0) != m) {
            stop("Number of rows of stat0 must equal length of stat.",
                call. = FALSE
            )
        }
        stat0 <- (stat0 - matrix(rep(stat, ncol(stat0)), byrow = FALSE, nrow = m)) >= 0
        p <- apply(stat0, 1, mean)
        p <- pmax(p, 1 / ncol(stat0))
    }
    return(p)
}
