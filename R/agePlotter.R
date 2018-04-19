#' Age plotter
#'
#' Plots a Y output by age
#'
#' @param y A vector with the Y outcome.
#' @param age A vector with the age information.
#' @param mod A model matrix.
#' @param mainText A character vector with the main title.
#' @param smoothIt A logical indicating whether to smooth Y.
#' @param jitter A logical indicating whether to add jitter on the X axis.
#' @param ageLabel The location of the age legend.
#' @param orderByAge A logical indicating whether to sort the observations by
#' age.
#' @param ylim A length two vector with the Y axis limits.
#' @param ageBreaks The age cutoffs for the different groups.
#' @param ylab The Y axis label.
#' @param pointColor An integer indicating the palette color to use for
#' the points.
#' @param lineColor An integer indicating the paletter color to use for the
#' lines.
#' @param alreadyFitted The output of \link[stats]{fitted} on a linear model
#' if you already calculated it. If so, `y` will be ignored.
#' @param ... Additional parameters to pass to \link[graphics]{plot}.
#' @details
#' \code{pointColor} can be a vector of length equal to \code{age} and have
#' multiple values in which case \code{lineColor} has to have a length
#' equal to the number of unique \code{pointColor} values. Specifying
#' this will draw a line for each unique \code{pointColor}.
#'
#' @return A nice plot =)
#'
#' @export
#' @author Andrew E Jaffe
#'
#' @import RColorBrewer
#' @import rafalib
#' @importFrom grDevices palette
#' @importFrom graphics axis layout lines mtext par plot text
#' @importFrom stats lm quantile
#'
#' @examples
#'

agePlotter <- function(y, age, mod = matrix(rep(1, length(y)), ncol=1),
    mainText, smoothIt = TRUE, jitter = TRUE, ageLabel = "bottom",
    orderByAge = TRUE, ylim = NULL, ageBreaks = c(-1, 0, 1, 10, 100),
    ylab = "Adjusted Expression", pointColor = 2, lineColor = 1,
    alreadyFitted = NULL, ...) {

    stopifnot(length(ageBreaks) >= 4)
    stopifnot(ageLabel %in% c('bottom', 'top'))
    stopifnot(length(lineColor) == length(unique(pointColor)))

    if(orderByAge) {
        oo <- order(age, decreasing = FALSE)
        y <- y[oo]
        age <- age[oo]
        mod <- mod[oo,, drop = FALSE]
        if(length(pointColor) == length(age)) pointColor <- pointColor[oo]
        if(!is.null(alreadyFitted)) alreadyFitted <- alreadyFitted[oo]
    }

    nfits <- length(unique(pointColor))
    if(is.null(alreadyFitted)) {
        if(nfits == 1) {
            fit <- fitted(lm(y ~ mod - 1))
        } else {
            fit <- lapply(unique(pointColor), function(col) {
                fitted(lm(y[pointColor == col] ~ mod[pointColor == col, , drop = FALSE] - 1))
            })
        }
    } else {
        if(nfits > 1) warning('Only one line will be draw. Try again specifying "mod"')
        nfits <- 1
        fit <- alreadyFitted
    }

    fetal <- cut(age, breaks = ageBreaks, lab = FALSE)
    fIndex <- splitit(fetal)

    make_line <- function(i, case0 = FALSE) {
        if(nfits == 1) {
            if(case0) {
                lines(age[fIndex[[i]]][-1], fit[fIndex[[i]]][-1], col = lineColor, lwd = 6)
            } else {
                lines(age[fIndex[[i]]], fit[fIndex[[i]]], col = lineColor, lwd = 6)
            }

        } else {
            for(j in seq_len(nfits)) {
                nfit_split <- splitit(cut(age[pointColor == unique(pointColor)[j]],
                    breaks = ageBreaks, lab = FALSE))
                if(case0) {
                    lines(age[pointColor == unique(pointColor)[j]][nfit_split[[i]]][-1],
                        fit[[j]][nfit_split[[i]]][-1],
                        col = lineColor[j], lwd = 6)
                } else {
                    lines(age[pointColor == unique(pointColor)[j]][nfit_split[[i]]], fit[[j]][nfit_split[[i]]],
                        col = lineColor[j], lwd = 6)
                }

            }
        }
    }

    nBreaks <- length(ageBreaks) - 1
    layout(matrix(rep(seq_len(nBreaks), c(5, rep(2, nBreaks - 2), 4)),
        nrow = 1, byrow = TRUE))
    palette(brewer.pal(8, "Set1"))

    par(mar = c(4, 5, 3, 0.45))
    if(is.null(ylim)) ylims <- range(y,na.rm=TRUE) else ylims <- ylim

    if(jitter) xx <- jitter(age, amount=0.005) else xx <- age
    plot(y ~ xx,
        subset = fIndex[[1]],
        main = "", ylab = ylab, xlab = "",
        ylim = ylims, cex.axis = 1.5, cex.lab = 1.75,
        pch = 21, cex = 1.4, bg = pointColor, xaxt = 'n',
        xlim = c(range(age[fIndex[[1]]]) + c(-0.01, 0.01)), ...)

    if(smoothIt) make_line(1)


    fetal_rang <- round(range(age[fIndex[[1]]]) * 52 + 40, 0)
    fetal_ax <- seq(fetal_rang[1], fetal_rang[2], round(diff(fetal_rang) / 4, 0))

    axis(1, at = (fetal_ax - 40) / 52,
        labels = fetal_ax, 1, cex.axis = 1.5)

    if(ageLabel == "bottom") {
        text(x = quantile(age[fIndex[[1]]], 0.33), y = min(ylims), "PCW",
            cex=1.5)
    } else if(ageLabel == "top") {
        text(x = quantile(age[fIndex[[1]]], 0.33), y = max(ylims), "PCW",
            cex=1.5)
    }

    # infant + child
    par(mar = c(4, 0.25,3,0.25))
    for(j in 2:(nBreaks - 1)) {
        plot(y ~ age, subset=fIndex[[j]],
            main = "", ylab = "" , xlab = "", yaxt = "n", cex = 1.4,
            xlim = range(age[fIndex[[j]]]) + c(-0.03, 0.03),
            ylim = ylims, cex.axis = 1.5, pch = 21, bg = pointColor, ...)

        if(ageBreaks[2] == 0 & smoothIt) make_line(j)
        if(ageBreaks[2] < 0 & smoothIt) make_line(j, case0 = TRUE)
    }

    # adults
    par(mar = c(4, 0.25,3,1))
    plot(y ~ age, subset=fIndex[[length(fIndex)]],
        main = "", ylab = "", xlab = "", yaxt = "n", cex = 1.4,
        xlim = range(age[fIndex[[length(fIndex)]]]) + c(-0.01, 0.01),
        ylim = ylims, cex.axis = 1.5, pch = 21, bg = pointColor, ...)

    if(smoothIt) make_line(length(fIndex))

    mtext(mainText, outer = TRUE, line = -2.5, cex = 1.35)
    mtext("Age", side = 1, outer = TRUE, line = -1.5, cex = 1.35)
}
