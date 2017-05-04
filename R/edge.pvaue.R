#' Function to compute empirical p-values 
#'
#' Given a pattern, it splits and returns the slot of interest.
#'
#' @param stat vector of observed statistics
#' @param stat0 matrix of null statistics, where rows are features and columns are iterations
#' @param pool features are interchangeable, e.g. unlist all null statistics
#' @param ... Additional arguments passed to \link[base]{strsplit}.
#'
#' @return A vector of empirical p-values.
#'
#' @export
#' @author Jeffrey T Leek
#'
#' @examples
#'
#'
edge.pvalue <- function(stat, stat0, pool=TRUE) {
  err.func <- "edge.pvalue"
  m <- length(stat)
  if(pool==TRUE) {
    if(is.matrix(stat0)) {stat0 <- as.vector(stat0)}
    m0 <- length(stat0) 
    v <- c(rep(T, m), rep(F, m0))
    v <- v[order(c(stat,stat0), decreasing = TRUE)]
    u <- 1:length(v)
    w <- 1:m
    p <- (u[v==TRUE]-w)/m0
    p <- p[rank(-stat)]
    p <- pmax(p,1/m0)
  } else {
    if(is.vector(stat0)) {
      err.msg(err.func,"stat0 must be a matrix.")
      return(invisible(1))
    }
    if(ncol(stat0)==m) {stat0 <- t(stat0)}
    if(nrow(stat0)!=m){
      err.msg(err.func,"Number of rows of stat0 must equal length of stat.")
      return(invisible(1))
    }
    stat0 <- (stat0 - matrix(rep(stat,ncol(stat0)),byrow=FALSE,nrow=m)) >= 0
    p <- apply(stat0,1,mean)
    p <- pmax(p,1/ncol(stat0))
  }
  return(p)
}