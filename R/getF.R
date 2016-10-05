#' Given two nested models, calculate the F-statistics
#'
#' For some data Y organized in a matrix, calculate a F-statitic per row
#' comparing two nested models.
#'
#' @param fit An object created with \link[limma]{lmFit} using the alternative
#' model (the larger model).
#' @param fit0 An object created with \link[limma]{lmFit} using the null
#' model (the smaller model, nested in the larger one).
#' @param theData A matrix with the data used in the \link[limma]{lmFit} call.
#' 
#'
#' @return A data.frame with the F-statistics (\code{fstat}), the degrees of
#' freedom for the nallternative model (\code{df1}), the null model
#' (\code{df0}), and the p-value given the F-distribution (\code{f_pval}).
#'
#' @export
#' @author Andrew E Jaffe
#' @import limma
#'
#' @examples
#'
#' set.seed(20161005)
#' ## From limma::limFit example page:
#' sd <- 0.3 * sqrt(4 / rchisq(100, df = 4))
#' y <- matrix(rnorm(100 * 6, sd = sd), 100, 6)
#' rownames(y) <- paste('Gene', seq_len(100))
#' y[1:2, 4:6] <- y[1:2, 4:6] + 4
#'
#' ## Define the alternative and null models
#' pheno <- data.frame(group = rep(c(0, 1), each = 3), RIN = runif(6) + 8)
#' mod <- model.matrix(~ pheno$group + pheno$RIN)
#' mod0 <- model.matrix(~ pheno$RIN)
#'
#' ## Fit the models
#' library('limma')
#' fit <- lmFit(y, mod)
#' fit0 <- lmFit(y, mod0)
#'
#' ## Calculate the F statistics for these nested models
#' finfo <- getF(fit, fit0, y)
#' head(finfo)
#'
#' ## You can then use p.adjust() for multiple testing corrections
#' qvals <- p.adjust(finfo$f_pval, 'fdr')
#' summary(qvals)
#'

getF <- function(fit, fit0, theData) {
    ## Check inputs
    stopifnot(is.matrix(theData))
    stopifnot(ncol(fit$design) > ncol(fit0$design))
    
    
	rss1 <- rowSums((fitted(fit) - theData)^2)
	df1 <- ncol(fit$coef) 
	rss0 <- rowSums((fitted(fit0) - theData)^2)
	df0 <- ncol(fit0$coef)
	fstat <- ((rss0 - rss1) / (df1 - df0)) / (rss1 / (ncol(theData) - df1))
	f_pval <- pf(fstat, df1 - df0, ncol(theData) - df1, lower.tail = FALSE)
	fout <- cbind(fstat, df1 - df0, ncol(theData) - df1, f_pval)
	colnames(fout)[2:3] <- c('df1', 'df0')
	fout <- data.frame(fout)
	return(fout)
}
