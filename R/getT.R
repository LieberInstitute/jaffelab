#' Get t-statistics from limma objects
#'
#' Get the t-statistics from a [lmFit][limma::lmFit] object.
#'
#' @param lmFitObj A [lmFit][limma::lmFit] object.
#' @param coiIndex An integer specifying the coefficent of interest in the model
#' for which to get the t-statistics for.
#'
#' @return A two column matrix with the t-statistics and the p-values of the
#' coefficient of interest.
#'
#' @export
#' @author Andrew E Jaffe
#' @importFrom methods is
#' @importFrom stats pt
#'
#' @examples
#'
#' ## From lmFit() example page
#' sd <- 0.3 * sqrt(4 / rchisq(100, df = 4))
#' y <- matrix(rnorm(100 * 6, sd = sd), 100, 6)
#' rownames(y) <- paste("Gene", 1:100)
#' y[1:2, 4:6] <- y[1:2, 4:6] + 2
#' design <- cbind(Grp1 = 1, Grp2vs1 = c(0, 0, 0, 1, 1, 1))
#'
#' # Ordinary fit
#' library("limma")
#' fit <- lmFit(y, design)
#' ## Get Ts
#' getT(fit)
getT <- function(lmFitObj, coiIndex = 2L) {
    ## Check inputs
    stopifnot(is(lmFitObj, "MArrayLM"))
    stopifnot(is.numeric(coiIndex) | is.integer(coiIndex))
    stopifnot(coiIndex >= 1 & coiIndex <= ncol(lmFitObj$coef))


    tt <- lmFitObj$coef[, coiIndex] / lmFitObj$stdev.unscaled[, coiIndex] /
        lmFitObj$sigma
    pp <- 2 * pt(-abs(tt), df = lmFitObj$df.residual)
    out <- cbind(tt, pp)
    colnames(out) <- c("Tstatistic", "pvalue")
    return(out)
}
