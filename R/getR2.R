#' Calculate R squared from a matrix
#'
#' Using limma, calculate the R squared from a matrix using nested models
#'
#' @param p A matrix that will be passed to [lmFit][limma::lmFit] `object`.
#' @param mod A model matrix for the alternative model (the larger one).
#' @param mod0 A model matrix for the null model (the smaller one). If `NULL`
#' then `p` will be used to calculate the residual sum of squares of the
#' null model.
#'
#' @return A data.frame with the R squared and the adjusted R squared.
#'
#' @export
#' @author Andrew E Jaffe, Leonardo Collado-Torres (examples)
#' @import limma
#'
#' @examples
#'
#' ## Define a model generating function for 30 'samples'
#' set.seed(20190827)
#' model_fun <- function(x) {
#'     ## Baseline + a group effect (2 groups) + a second covariate effect
#'     rnorm(30) +
#'         c(rnorm(15, mean = 3), rnorm(15, mean = 1)) +
#'         c(
#'             rnorm(5, sd = 0.5), rnorm(5, sd = 0.2, mean = 0.5),
#'             rnorm(5, sd = 0.2, mean = 0.9)
#'         )
#' }
#'
#' ## Generate the data for 20 'genes'
#' p <- t(sapply(1:20, model_fun))
#'
#' ## Define the phenotype data for these 30 'samples'
#' pheno <- data.frame(
#'     group = rep(c("A", "B"), each = 15),
#'     batch = rep(1:3, each = 5)
#' )
#'
#' ## Define a full model
#' mod <- with(pheno, model.matrix(~ group + batch))
#' ## and compute the R2 for each 'gene'
#' getR2(p, mod)
#'
#' ## Define a smaller model
#' mod0 <- with(pheno, model.matrix(~group))
#' ## And now compute the new R2 for each 'gene'
#' getR2(p, mod, mod0)
getR2 <- function(p, mod, mod0 = NULL) {
    fit1 <- lmFit(p, mod)
    rss1 <- rowSums((p - fitted(fit1))^2)
    n <- ncol(p)
    k <- ncol(mod) - 1

    if (is.null(mod0)) {
        rss0 <- rowSums((p - rowMeans(p))^2)
    } else {
        fit0 <- lmFit(p, mod0)
        rss0 <- rowSums((p - fitted(fit0))^2)
    }

    r2 <- 1 - (rss1 / rss0)
    r2adj <- 1 - ((1 - r2) * (n - 1)) / (n - k - 1)
    out <- data.frame(R2 = r2, Adjusted_R2 = r2adj)
    return(out)
}
