#' Given a set of regions and bigwig files, compute the coverage matrix for
#' using \code{bwtool} and build a
#' \link[SummarizedExperiment]{RangedSummarizedExperiment-class} object.
#'
#' Given a set of regions and bigwig files, compute the coverage matrix for
#' using \code{bwtool} and build a
#' \link[SummarizedExperiment]{RangedSummarizedExperiment-class} object.
#'
#' @param bws A named vector with the paths to the bigWig files. The names
#' are used as sample ids.
#' @param regions A \link[GenomicRanges]{GRanges-class} object with regions
#' for which to calculate the coverage matrix.
#' @param strand Either *, + or -. If set to * (default) then all regions are
#' used. Otherwise the matrix is subset to the regions of the corresponding
#' strand. The users should supply the correct corresponding list of bigWig
#' files in \code{bws}.
#' @param pheno \code{NULL} by default. Specify the data.frame with the same
#' length as \code{bws} to be used in the resulting RSE object.
#' @param bwtool The path to \code{bwtool}. Uses as the default the
#' location at JHPCE.
#' @param bpparam A \link[BiocParallel]{BiocParallelParam-class} instance which
#' will be used to calculate the coverage matrix in parallel. By default, 
#' \link[BiocParallel]{SerialParam-class} will be used.
#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
#' @param sumsdir The path to an existing directory where the \code{bwtool}
#' sum tsv files will be saved. We recommend setting this to a value beyond
#' the default one. Use separate output directories per strand.
#' @param commands_only If \code{TRUE} the bwtool commands will be saved in a
#' file called coverage_bwtool_strandSTRAND.txt and exit without running
#' \code{bwtool}. This is useful if you have a very large regions set and want
#' to run the commands in an array job. Then run
#' \code{coverage_matrix_bwtool(commands_only = FALSE)} to create the RSE
#' object(s).
#' 
#'
#' @return A \link[SummarizedExperiment]{RangedSummarizedExperiment-class}
#' object with the counts stored in the assays slot. 
#'
#' @details Based on \link[recount.bwtool]{coverage_matrix_bwtool}.
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @importFrom methods is
#' @import GenomicRanges BiocParallel rtracklayer SummarizedExperiment S4Vectors
#' recount.bwtool
#'
#' @seealso \link[recount.bwtool]{coverage_matrix_bwtool}
#'
#' @examples
#' if(.Platform$OS.type != 'windows') {
#' ## Disable the example for now. I'd have to figure out how to install
#' ## bwtool on travis
#' if(FALSE) {
#'     ## Reading BigWig files is not supported by rtracklayer on Windows
#'     ## (only needed for defining the regions in this example)
#'     
#'     ## TODO
#' }
#' }
#'

coverage_bwtool <- function(bws, regions, strand = '*', pheno = NULL,
    bwtool = '/dcl01/leek/data/bwtool/bwtool-1.0/bwtool',
    bpparam = NULL, verbose = TRUE, sumsdir = tempdir(),
    commands_only = FALSE) {
        
    ## Check inputs
    stopifnot(!is.null(names(bws)))
    stopifnot(strand %in% c('*', '+', '-'))
    stopifnot(all(file.exists(bws)))
    
    ## Build a pheno if missing
    if(is.null(pheno)) {
        pheno <- data.frame(bigwig_path = bws, bigwig_file = basename(bws),
            sample = names(bws))
    } else {
        stopifnot(nrow(pheno) != length(bws))
    }
    dir.create(sumsdir, recursive = TRUE, showWarnings = FALSE)
    
    ##  Subset regions by strand if needed
    if(strand != '*') {
        regions <- regions[strand(regions) == strand]
    }
    
    ## Export regions to a BED file if necessary
    bed <- file.path(sumsdir, paste0('coverage_bwtool_strand', strand, '-',
        Sys.Date(), '.bed'))
    
    if(!file.exists(bed)) {
        if (verbose) message(paste(Sys.time(), 'creating the BED file', bed))
        rtracklayer::export(regions, con = bed, format='BED')
        stopifnot(file.exists(bed))
    }
    
    ## Define bpparam
    if(is.null(bpparam)) bpparam <- BiocParallel::SerialParam()
    
    ## Run bwtool and load the data
    counts <- bpmapply(recount.bwtool:::.run_bwtool, bws, names(bws),
        MoreArgs = list('bwtool' = bwtool, 'bed' = bed, 'sumsdir' = sumsdir,
        'verbose' = verbose, 'commands_only' = commands_only),
        SIMPLIFY = FALSE, BPPARAM = bpparam)
    if(commands_only) {
        commands <- unlist(counts)
        cat(commands, file = paste0('coverage_bwtool_strand', strand,
            '.txt'), sep = '\n')
        return(invisible(NULL))
    }
    
    ## Group results from all files
    counts <- do.call(cbind, counts)

    ## Build a RSE object
    rse <- SummarizedExperiment(assays = list('counts' = counts),
            colData = DataFrame(pheno), rowRanges = regions)
            
    ## Finish
    return(rse)
}
