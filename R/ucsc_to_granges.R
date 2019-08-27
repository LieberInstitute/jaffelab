#' Convert UCSC coordinates to a GRanges object
#'
#' This function takes a string of UCSC genome browser coordinates
#' and returns them as a GRanges object.
#'
#' @param coords A character vector of UCSC genome browser coordinates.
#'
#' @return A GRanges object.
#' @seealso [granges_to_ucsc]
#' @author Leonardo Collado-Torres, Andrew E Jaffe
#' @export
#' @importFrom GenomicRanges GRanges
#'
#' @examples
#' ucsc_coords <- c('chr1:1000-2000')
#' ucsc_to_granges(ucsc_coords)
#'
#'
ucsc_to_granges <- function(coords) {
    stopifnot(is.character(coords))
    chr <- ss(coords, ":")
    start <- as.numeric(ss(ss(coords, ":", 2), "-"))
    end <- as.numeric(ss(ss(coords, ":", 2), "-", 2))

    GRanges(chr, IRanges(start, end))
}


#' Convert a GRanges object to UCSC coordinates
#'
#' This function takes a GRanges objects and returns a string of UCSC genome
#' browser coordinates.
#'
#' @param gr A GRanges object.
#'
#' @return A character vector of UCSC genome browser coordinates.
#' @seealso [ucsc_to_granges]
#' @author Leonardo Collado-Torres
#' @export
#' @importFrom GenomicRanges GRanges
#'
#' @examples
#' ucsc_coords <- c('chr1:1000-2000')
#' granges_to_ucsc(ucsc_to_granges(ucsc_coords))
#'
#'
granges_to_ucsc <- function(gr) {
    stopifnot(is(gr, 'GRanges'))

    chr <- as.character(seqnames(gr))
    start <- start(gr)
    end <- end(gr)

    paste0(chr, ':', start, '-', end)
}
