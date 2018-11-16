#' Count junctions
#'
#' Count junctions from either TopHat2, HISAT2 or Rail-RNA output
#'
#' @param junctionFiles A character vector with the full paths to the junction
#' files. Can alternatively be a list of \link[GenomicRanges]{GRanges-class}
#' objects with the junction count information.
#' @param sampleNames A character vector of the same length as `junctionsFiles`
#' to use as the sample names.
#' @param output Either `Count` (TopHat2, HISAT2) or `Rail` (Rail-RNA).
#' @param minOverhang minimum overhang.
#' @param strandSpecific a logical specifying whether the library is strand
#' specific or not.
#' @param illuminaStranded a logical indicating whether the Illumina library
#' is stranded or not.
#' @param minCount Minimum count.
#' @param maxCores The maximum number of cores to use. By default one.
#'
#' @return A two element list with a `DataFrame` and a
#' \link[GenomicRanges]{GRanges-class} object with the counts and the
#' annotation used.
#'
#' @export
#' @author Andrew E Jaffe
#' @import GenomicRanges
#' @import parallel
#' @importFrom utils read.delim
#' @importFrom IRanges IRanges
#' @import S4Vectors
#'
#' @examples
#'

junctionCount<- function(junctionFiles, sampleNames = names(junctionFiles),
    output = c("Count", "Rail"), minOverhang = 0,
    strandSpecific = FALSE, illuminaStranded = FALSE,
    minCount = 1, maxCores = 1) {

    stopifnot(length(junctionFiles) == length(sampleNames))
    stopifnot(output %in% c('Count', 'Rail'))

    names(junctionFiles) <- sampleNames
    message(paste(Sys.time(), 'reading in data'))

    if(all(is.character(junctionFiles))) {
        theData <- mclapply(junctionFiles, function(x) {
            if(output == "Rail") {
                y <- read.delim(x, skip = 1, header = FALSE,
                    colClasses = c("character", "integer",
                    "integer", "integer", "integer", "integer"),
                    col.names = c("chr", "start", "end", "leftHang",
                    "rightHang", "count"))
                y <- y[y$count >= minCount,] # filter based on min number
                y <- y[y$leftHang > minOverhang & y$rightHang > minOverhang,]
            } else if(output == "Count") {
                y <- read.delim(x, skip = 1, header = FALSE,
                col.names = c("chr", "start","end", "strand", "count"),
                colClasses = c("character", "integer", "integer",
                "character","integer"))
                y <- y[y$count >= minCount, ] # filter based on min number
                weird <- which(y$strand=="?")
                if(length(weird) > 0) y <- y[-weird, ]
            }

            gr <- GRanges(y$chr, IRanges(y$start, y$end), strand = y$strand,
                count = y$count)
            return(gr)
        }, mc.cores = maxCores)
    } else {
        theData <- junctionFiles
        stopifnot(all(sapply(theData, class) == "GRanges"))
    }
    message(paste(Sys.time(), 'creating master table of junctions'))

    ## turn into GRangesList
    ### THIS STEP IS SLOW...
    grList <- GRangesList(theData)

    # each dataset should be checked
    if(illuminaStranded & strandSpecific) {
        grList <- GRangesList(mclapply(grList, function(x) {
            strand(x) = ifelse(strand(x) == "+", "-","+")
            return(x)
        }, mc.cores = maxCores))
    }

    ## get into GRanges object of unique junctions
    fullGR <- unlist(grList)
    if(!strandSpecific) strand(fullGR) <- "*"

    fullGR <- fullGR[!duplicated(fullGR)] # or unique(fullGR)
    fullGR <- sort(fullGR)
    fullGR$count <- NULL

    message(paste(Sys.time(), 'there are', length(fullGR), 'total junctions'))
    message(paste(Sys.time(), 'populating count matrix'))

    jNames <- paste0(as.character(seqnames(fullGR)), ":", start(fullGR), "-",
        end(fullGR), "(", as.character(strand(fullGR)), ")")

    ## match GRanges
    options(warn = -1)
    mList <- mclapply(grList, match, fullGR,
        ignore.strand = !strandSpecific, mc.cores = maxCores)
    options(warn = 0)

    countList <- mList # initiate
    M <- length(jNames)

    ## fill in matrix
    message(paste(Sys.time(), 'filling in the count matrix'))
    for(i in seq(along=grList)) {
        if(i %% 25 == 0) cat(".")
        cc <- rep(0,M)
        cc[mList[[i]]] <- theData[[i]]$count
        countList[[i]] <- Rle(cc)
    }
    countDF <- DataFrame(countList, row.names = jNames, check.names = FALSE)

    names(fullGR) <- jNames
    ## return matrix and GRanges object
    out <- list(countDF = countDF, anno = fullGR)
    return(out)
}
