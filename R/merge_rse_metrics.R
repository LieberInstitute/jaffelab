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
#' 'mitoMapped', 'totalMapped', 'ERCCsumLogErr' as either NumericList() or
#' IntegerList() objects in the colData() columns.
#'
#' @return A [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#' with merged columns for 'concordMapRate', 'overallMapRate', 'mitoRate',
#' 'rRNA_rate', 'totalAssignedGene', 'numMapped', 'numReads', 'numUnmapped',
#' 'mitoMapped', 'totalMapped', 'ERCCsumLogErr'
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
#'     feature_id = sprintf("ID%03d", seq_len(200))
#' )
#'
#' ## Function for generating some data
#' rand_gen <- function(mu = 1e6, sd = 1e3, digits = 0) {
#'     IRanges::NumericList(lapply(seq_len(6), function(x) {
#'         round(rnorm(2, mu, sd), digits)
#'     }))
#' }
#' colData <- DataFrame(
#'     Treatment = rep(c("ChIP", "Input"), 3),
#'     row.names = LETTERS[seq_len(6)],
#'     concordMapRate = rand_gen(0.5, 0.05, 10),
#'     overallMapRate = rand_gen(0.7, 0.05, 10),
#'     mitoRate = rand_gen(0.2, 0.05, 10),
#'     rRNA_rate = rand_gen(0.1, 0.01, 10),
#'     totalAssignedGene = rand_gen(0.5, 0.05, 10),
#'     numMapped = rand_gen(),
#'     numReads = rand_gen(1e6 + 6e5),
#'     numUnmapped = rand_gen(5e5),
#'     mitoMapped = rand_gen(1e5),
#'     totalMapped = rand_gen(1e6 + 1e5),
#'     ERCCsumLogErr = rand_gen(-28.5, 4, 6)
#' )
#' rse <- SummarizedExperiment(
#'     assays = SimpleList(counts = counts),
#'     rowRanges = rowRanges, colData = colData
#' )
#' colData(rse)
#' rse_merged <- merge_rse_metrics(rse)
#' colData(rse_merged)
#'
#' ## Use the BSP2 real data
#' local_bsp2 <- tempfile("BSP2_rse_gene.Rdata")
#' download.file(
#'     "https://s3.us-east-2.amazonaws.com/libd-brainseq2/rse_gene_unfiltered.Rdata",
#'     destfile = local_bsp2,
#'     mode = "wb"
#' )
#' load(local_bsp2, verbose = TRUE)
#' # lobstr::obj_size(rse_gene)
#' # 872.12 MB
#' bsp2_rse_merged <- merge_rse_metrics(rse_gene)
#' colData(bsp2_rse_merged)
#'
#' ## Compare the ERCCsumLogErr approximation against the mean of the
#' ## error values. Note that you would actually have to load the actual
#' ## ERCC quantification output files in order to compute the correct
#' ## ERCCsumLogErr values.
#' bsp2_rse_merged$mean_ERCCsumLogErr <-
#'     sapply(rse_gene$ERCCsumLogErr, mean)
#' with(
#'     colData(bsp2_rse_merged)[lengths(rse_gene$ERCCsumLogErr) > 1, ],
#'     plot(
#'         mean_ERCCsumLogErr,
#'         ERCCsumLogErr,
#'         ylab = "Approximate ERCCsumLogErr"
#'     )
#' )
#' abline(a = 0, b = 1, col = "red")
#' with(
#'     colData(bsp2_rse_merged)[lengths(rse_gene$ERCCsumLogErr) > 1, ],
#'     summary(mean_ERCCsumLogErr - ERCCsumLogErr)
#' )
#'
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

    if ("ERCCsumLogErr" %in% colnames(colData(rse))) {
        message(Sys.time(), " attempting to approximate ERCCsumLogErr.")
        stopifnot(is(rse$ERCCsumLogErr, "NumericList"))

        # ## Read in the ERCC expected concentration
        # spikeIns <- read.delim("https://raw.githubusercontent.com/LieberInstitute/SPEAQeasy/master/Annotation/ERCC/ercc_actual_conc.txt",
        #     as.is = TRUE, row.names = 2
        # )
        #
        # ## Identify the number of ERCC sequences
        # n_erccs <- nrow(spikeIns)
        #
        # ## Calculate the ERCC expected concentration
        # ercc_expected_conc <- sum(log2(10 * spikeIns[, "concentration.in.Mix.1..attomoles.ul."] + 1))
        # dput(ercc_expected_conc)
        # dput(n_erccs)
        #
        ## Use pre-computed values
        ercc_expected_conc <- 643.072778500804
        n_erccs <- 92L

        ## Then compute the approximate result starting from a NumericList input
        rse$ERCCsumLogErr <- sapply(rse$ERCCsumLogErr, function(sumlogerr) {
            log2(mean((2^(
                sumlogerr + ercc_expected_conc
            ))^(1 / n_erccs))^n_erccs) - ercc_expected_conc
        })
    }

    return(rse)
}
