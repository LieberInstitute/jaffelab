#' corner, an improved head()
#'
#' Returns the first part of an object similar to head() but subsets
#' columns as well as rows, and also handles lists of objects.
#'
#' @param mat A DataFrame, data.frame, GRanges, matrix, tibble, or list
#' @param n Integer number of rows and columns you want returned
#'
#' @return An n x n object of the original class (or a list of n nxn objects)
#' @export
#' @author Emily E. Burke
#'
#' @examples
#'
#' m <- matrix(rnorm(1000), ncol = 10)
#' corner(m)
#'
#' lis <- list(iris, mtcars, matrix(rnorm(1000), ncol = 10))
#' corner(lis)
#'
corner = function(mat, n=6) {
	types = c("list","DataFrame","data.frame","GRanges","matrix","tbl_df")
	if(!(class(mat)[1] %in% types)) {
		stop(paste0('The class of your object (',class(mat),') is not handled by this function.')) }
	if (class(mat)[1] == "list") {
		# return corners of first n items of list
		return(lapply(mat[seq_len(min(length(mat), n))], function(x) x[seq_len(min(nrow(x), n)),seq_len(min(ncol(x), n))]) )
	} else {
		return(mat[seq_len(min(nrow(mat), n)),seq_len(min(ncol(mat), n))])
	}
}


