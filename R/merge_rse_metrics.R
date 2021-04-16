#' Merge sample metadata across several columns
#'
#' This function merges the sample metadata from a RangedSummarizedExperiment
#' object for some of the variables that we (at Jaffelab) typically create
#' with our RNA-seq processing pipeline for samples that are sequenced
#' in more than one lane.
#'
#' @param rse A [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#' object with 'concordMapRate', 'overallMapRate', 'mitoRate',
#' 'rRNA_rate', 'totalAssignedGene', 'numMapped', 'numReads', 'numUnmapped',
#' 'mitoMapped', 'totalMapped' as either NumericList() or IntegerList() objects
#' in the colData() columns.
#'
#' @return A [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#' with merged columns for 'concordMapRate', 'overallMapRate', 'mitoRate',
#' 'rRNA_rate', 'totalAssignedGene', 'numMapped', 'numReads', 'numUnmapped',
#' 'mitoMapped', 'totalMapped'.
#'
#' @export
#' @author Leonardo Collado-Torres, Andrew E Jaffe
#' @import SummarizedExperiment
#'
#' @examples
#'
#' ## Taken from ?SummarizedExperiment
#' library("SummarizedExperiment")
#' nrows <- 200
#' ncols <- 6
#' set.seed(20181116)
#' counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#' rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
#'     IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
#'     strand = sample(c("+", "-"), 200, TRUE),
#'     feature_id = sprintf("ID%03d", 1:200)
#' )
#'
#' ## Function for generating some data
#' rand_gen <- function(mu = 1e6, sd = 1e3, digits = 0) {
#'     IRanges::NumericList(lapply(1:6, function(x) {
#'         round(rnorm(2, mu, sd), digits)
#'     }))
#' }
#' colData <- DataFrame(
#'     Treatment = rep(c("ChIP", "Input"), 3),
#'     row.names = LETTERS[1:6],
#'     concordMapRate = rand_gen(0.5, 0.05, 10),
#'     overallMapRate = rand_gen(0.7, 0.05, 10),
#'     mitoRate = rand_gen(0.2, 0.05, 10),
#'     rRNA_rate = rand_gen(0.1, 0.01, 10),
#'     totalAssignedGene = rand_gen(0.5, 0.05, 10),
#'     numMapped = rand_gen(),
#'     numReads = rand_gen(1e6 + 6e5),
#'     numUnmapped = rand_gen(5e5),
#'     mitoMapped = rand_gen(1e5),
#'     totalMapped = rand_gen(1e6 + 1e5)
#' )
#' rse <- SummarizedExperiment(
#'     assays = SimpleList(counts = counts),
#'     rowRanges = rowRanges, colData = colData
#' )
#' colData(rse)
#' rse_merged <- merge_rse_metrics(rse)
#' colData(rse_merged)
merge_rse_metrics <- function(rse) {
    stopifnot(is(rse, "RangedSummarizedExperiment"))
    stopifnot(
        c(
            "concordMapRate", "overallMapRate", "mitoRate", "rRNA_rate",
            "totalAssignedGene", "numMapped", "numReads", "numUnmapped",
            "mitoMapped", "totalMapped"
        ) %in%
            colnames(SummarizedExperiment::colData(rse))
    )

    stopifnot(all(sapply(c(
        "concordMapRate", "overallMapRate", "mitoRate", "rRNA_rate",
        "totalAssignedGene", "numMapped", "numReads", "numUnmapped",
        "mitoMapped", "totalMapped"
    ), function(var) {
        is(colData(rse)[, var], "List")
    })))

    rse$concordMapRate <- mapply(function(r, n) {
        sum(r * n) / sum(n)
    }, rse$concordMapRate, rse$numReads)
    rse$overallMapRate <- mapply(function(r, n) {
        sum(r * n) / sum(n)
    }, rse$overallMapRate, rse$numReads)
    rse$mitoRate <- mapply(function(r, n) {
        sum(r * n) / sum(n)
    }, rse$mitoRate, rse$numMapped)
    rse$rRNA_rate <- mapply(function(r, n) {
        sum(r * n) / sum(n)
    }, rse$rRNA_rate, rse$numMapped)
    rse$totalAssignedGene <- mapply(function(r, n) {
        sum(r * n) / sum(n)
    }, rse$totalAssignedGene, rse$numMapped)

    rse$numMapped <- sapply(rse$numMapped, sum)
    rse$numReads <- sapply(rse$numReads, sum)
    rse$numUnmapped <- sapply(rse$numUnmapped, sum)
    rse$mitoMapped <- sapply(rse$mitoMapped, sum)
    rse$totalMapped <- sapply(rse$totalMapped, sum)
    return(rse)
}
