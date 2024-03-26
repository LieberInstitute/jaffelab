#' Find expression cutoffs
#'
#' This function finds expression cutoffs based on estimating the second
#' breakpoint using [segmented][segmented::segmented] applied to two curves. The
#' first curve is the number of features passing the cutoff assuming an initial
#' large drop, an intermediate section with moderate drop, and then a stable
#' tail. The second curve is the mean number of expressed samples (non-zero
#' expression) for all genes at each given cutoff. This curve increases rapidly
#' then has a section with a moderate increase, and finally a long tail. Using
#' [segmented][segmented] we find that second breakpoint. The suggested
#' expression cutoff returned is the average of the two suggested cutoffs.
#'
#' @param expr A matrix with the expression values with features in the rows and
#' samples in the columns in RPKM format (typically). `expr` can be any matrix
#' as `expression_cutoff()` will work with the row means of the matrix. You'll
#' most likely want to use a normalized expression matrix, although this
#' function will work with any matrix.
#' @param max_cut Maximum expression cutoff to consider. Increasing this value
#' does then to increase the suggested cutoffs. If set to `NULL`, this
#' will be chosen automatically from 1 to 5.
#' @param seed Set this argument for increased reproducibility of the results.
#' This is passed to [seg.control][segmented::seg.control]. We highly recommend
#' specifiying this value.
#' @param n.boot This argument controls the stability of the suggested cutoffs.
#' It is passed to [seg.control][segmented::seg.control].
#' @param k The number of breakpoints to consider. If you increase this value
#' the estimated breakpoints are less stable.
#'
#'
#' @return A two element vector with the suggested cutoffs based on the
#' number of features expressed and the mean number of expressed samples
#' across all features. The function also makes plots that visualize
#' both curves and show the identified breakpoints as vertical grey lines.
#'
#' @export
#' @author Leonardo Collado-Torres
#' @import segmented
#' @importFrom graphics abline boxplot legend
#' @importFrom stats median
#'
#' @examples
#'
#' ## Simulate expression data. Based on the limma::lmFit() man page
#' sd <- 0.3 * sqrt(4 / rchisq(100, df = 4))
#' y <- matrix(rnorm(100 * 6, sd = sd), 100, 6)
#' rownames(y) <- paste("Gene", 1:100)
#' y[1:2, 4:6] <- y[1:2, 4:6] + 2
#'
#' ## Make the expression data look like RPKM values
#' y2 <- apply(y, 2, function(z) {
#'     res <- z + abs(min(z))
#'     res[sample(seq_len(100), 10)] <- 0
#'     res
#' })
#' summary(y2)
#'
#' expression_cutoff(y2)
#'
#' ## Or same the plots to a PDF file
#'
#' pdf("test_expression_cutoff.pdf", width = 12)
#' expression_cutoff(y2)
#' dev.off()
#'
#' ## View the pdf with the following code
#' utils::browseURL("test_expression_cutoff.pdf")
expression_cutoff <- function(expr, max_cut = 1, seed = NULL, n.boot = 2000,
    k = 2) {
    meanExpr <- rowMeans(expr)


    if (is.null(max_cut)) {
        cuts <- seq_len(500) / 100



        cuts_m <- which(sapply(cuts, function(cut) {
            mean(meanExpr > cut) <= 0.25
        }))[1]
        cuts <- c(-0.01, 0, seq_len(cuts_m) / 100)
    } else {
        cuts <- c(-0.01, 0, seq_len(round(max_cut * 100)) / 100)
    }


    cuts_df <- data.frame(do.call(rbind, lapply(cuts, function(cut) {
        c(
            "cut" = cut,
            "features_n" = sum(meanExpr > cut),
            "features_percent" = mean(meanExpr > cut)
        )
    })))

    cutl <- seq_len(length(cuts))
    f <- lm(features_percent ~ cutl, data = cuts_df)
    t.l <- cutl[length(cutl)]

    seg <- segmented(f,
        seg.Z = ~cutl,
        psi = round(seq(1, t.l, length.out = k + 2)[2:(k + 1)]),
        control = seg.control(n.boot = n.boot, seed = seed)
    )
    # round(seg$psi[, 2])
    # cuts[round(seg$psi[, 2])]
    # plot(seg)

    plot(cuts_df$cut, cuts_df$features_n,
        type = "o", col = "purple",
        ylab = "Number of features > cutoff", xlab = "mean expression cutoff",
        pch = 20
    )
    abline(v = cuts[round(seg$psi[, 2])], col = "grey80")


    plot(cuts_df$cut, cuts_df$features_percent * 100,
        type = "o",
        col = "purple", ylab = "Percent of features > cutoff",
        xlab = "mean expression cutoff", pch = 20
    )
    abline(v = cuts[round(seg$psi[, 2])], col = "grey80")
    abline(h = 25, col = "grey80")


    nonzero <- rowSums(expr > 0)

    nonzero_samples <- lapply(cuts, function(cut) {
        nonzero[meanExpr > cut]
    })
    names(nonzero_samples) <- cuts



    df2 <- data.frame(nonzeromean = sapply(nonzero_samples, mean))
    f2 <- lm(nonzeromean ~ cutl, data = df2)
    seg2 <- segmented(f2,
        seg.Z = ~cutl,
        psi = round(seq(1, t.l, length.out = k + 2)[2:(k + 1)]),
        control = seg.control(n.boot = n.boot, seed = seed)
    )
    # plot(seg2)
    # round(seg2$psi[, 2])
    # cuts[round(seg2$psi[, 2])]

    boxplot(nonzero_samples,
        xlab = "mean expression cutoff",
        ylab = "Number of expressed samples among features > cutoff",
        outline = FALSE
    )

    lines(seq_len(length(cuts)), sapply(nonzero_samples, mean),
        type = "o",
        col = "purple"
    )
    lines(seq_len(length(cuts)), sapply(nonzero_samples, median),
        type = "l",
        col = "orange"
    )
    abline(v = round(seg2$psi[, 2]), col = "grey80")
    legend("bottomright", c("mean", "median"),
        col = c("purple", "orange"),
        lwd = 2, bty = "n"
    )

    suggested <- c(
        "percent_features_cut" = cuts[round(seg$psi[, 2])][2],
        "samples_nonzero_cut" = cuts[round(seg2$psi[, 2])][2]
    )

    message(paste(
        Sys.time(), "the suggested expression cutoff is",
        round(mean(suggested), 2)
    ))
    return(suggested)
}
