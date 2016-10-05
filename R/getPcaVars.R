#' Calculates the percent of variance explained for principal components
#'
#' Given a PCA object created with \link[stats]{prcomp}, this function computes 
#' the percent of variance explained by each of the principal components.
#'
#' @param pca An object created with \link[stats]{prcomp}.
#' @param digits The number of significant digits to round to.
#'
#' @return A vector with the percent of variance explained for each of the
#' principal components in decreasing order.
#'
#' @export
#' @importFrom methods is
#' @author Andrew E Jaffe
#'
#' @examples
#'
#' pca <- prcomp(USArrests)
#' getPcaVars(pca)
#'

getPcaVars <- function(pca, digits = 3) {
    stopifnot(is(pca) == 'prcomp')
    signif(pca$sdev^2 / sum(pca$sdev^2) * 100, digits = digits)
}
