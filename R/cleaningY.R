#' Regress out covariates
#'
#' Regress out covariates such as surrogate variables or principal components.
#'
#' @param y A matrix such as the outcome matrix from sva or a gene expression
#' matrix.
#' @param mod A full rank model matrix.
#' @param P The number of SVs or PCs to protect based on the column order.
#' For example, `P=2` would keep the intercept term and a case vs diagnosis
#' term in a model that is ~ Dx + more covariates.
#'
#' @return An object of the same type as `y` with the SVs/PCs regressed out.
#'
#' @export
#' @author Rafael Irizarry, Leonardo Collado-Torres (examples)
#'
#' @examples
#'
#' ## Define a model generating function for 30 'samples'
#' set.seed(20190827)
#' model_fun <- function(x) {
#'     ## Baseline + a group effect (2 groups) + a second covariate effect
#'     rnorm(30) +
#'     c(rnorm(15, mean = 3), rnorm(15, mean = 1)) +
#'     c(rnorm(5, sd = 0.5), rnorm(5, sd = 0.2, mean = 0.5),
#'         rnorm(5, sd = 0.2, mean = 0.9))
#' }
#'
#' ## Generate the data for 20 'genes'
#' y <- t(sapply(1:20, model_fun))
#'
#' ## Define the phenotype data for these 30 'samples'
#' pheno <- data.frame(
#'     group = rep(c('A', 'B'), each = 15),
#'     batch = rep(1:3, each = 5)
#' )
#'
#' ## Define a full model
#' mod <- with(pheno, model.matrix(~ group + batch))
#'
#' ## Check the raw data for gene 1
#' boxplot(y[1, ] ~ pheno$group, ylab = 'Gene 1 Raw Expr')
#'
#' ## Now regress out the batch covariate from the gene expression matrix
#' y_clean_p2 <- cleaningY(y, mod, P = 2)
#'
#' ## Check the cleaned data for gene 1 (with P = 2)
#' boxplot(y_clean_p2[1, ] ~ pheno$group, ylab = 'Gene 1 Clean Expr (P = 2)')
#'
#' ## Or regress out the group and batch effects
#' y_clean_p3 <- cleaningY(y, mod, P = 1)
#'
#' ## Check the cleaned data for gene 1 (with P = 3)
#' boxplot(y_clean_p3[1, ] ~ pheno$group, ylab = 'Gene 1 Clean Expr (P = 3)')
#'
#'
#'
#'
#' ## The function also supports NAs observations as detailed below
#'
#' ## Make one observation 0, clean the data
#' y[1, 1] <- 0
#' y_clean_p2_0 <- cleaningY(y, mod, P = 2)
#' ## then NA and clean again
#' y[1, 1] <- NA
#' y_clean_p2_NA <- cleaningY(y, mod, P = 2)
#'
#' ## Compare the results
#' corner(y_clean_p2_0)
#' corner(y_clean_p2_NA)
#'
#' ## They are identical except for that NA in [1, 1]
#' table(y_clean_p2_0 -  y_clean_p2_NA, useNA = 'ifany')
#'
#' ## Compared to the original y, there are differences since we lost
#' ## one observation which affects all of the first row of the cleaned Y
#' y_clean_p2[1, ] - y_clean_p2_NA[1, ]
#' all(y_clean_p2[-1, ] - y_clean_p2_NA[-1, ] == 0)
#'

cleaningY <- function(y, mod, P) {
    stopifnot(P <= ncol(mod))
    Hat <- solve(t(mod) %*% mod) %*% t(mod)
    ## For dealing with NAs
    ## https://stackoverflow.com/questions/16535084/matrix-multiplication-with-scattered-na-values
    ty <- t(y)
    ty[is.na(ty)] <- 0
    beta <- (Hat %*% ty)
    ## Note that y might still have the NAs, and NA - a number = NA
    ## so there's no need to reset the NAs back on cleany
    cleany <- y - t(as.matrix(mod[, -c(seq_len(P))]) %*% beta[-seq_len(P), ])
    return(cleany)
}
