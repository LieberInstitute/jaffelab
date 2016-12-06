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
#'
#' @return A nice plot =)
#'
#' @export
#' @author Andrew E Jaffe
#'
#' @import RColorBrewer
#' @import rafalib
#'
#' @examples
#'

agePlotter <- function(y, age, mod = matrix(rep(1, length(y)), ncol=1),
    mainText, smoothIt = TRUE, jitter = TRUE, ageLabel = "bottom",
    orderByAge = TRUE, ylim = NULL, ageBreaks = c(-1, 0, 1, 10, 100), 
    ylab = "Adjusted Expression", pointColor = 2, lineColor = 1,
    alreadyFitted = NULL, ...) {
    
    if(orderByAge) {
        oo <- order(age, decreasing = FALSE)
        y <- y[oo]
        age <- age[oo]
        mod <- mod[oo,]
        if(!is.null(alreadyFitted)) alreadyFitted <- alreadyFitted[oo]
    }
    
    if(is.null(alreadyFitted)) {
        fit <- fitted(lm(y ~ mod - 1))
    } else fit <- alreadyFitted
    
    fetal <- cut(age, breaks = ageBreaks, lab = FALSE)
    fIndex <- splitit(fetal)    
    
    layout(matrix(rep(1:4, c(5, 2, 2, 4)), nr = 1, byrow = TRUE))
    palette(brewer.pal(8, "Set1"))
    
    par(mar = c(4, 5, 3, 0.45))
    if(is.null(ylim)) ylims <- range(y,na.rm=TRUE) else ylims <- ylim

    if(jitter) xx <- jitter(age, amount=0.005) else xx <- age
    plot(y ~ xx,
        subset = fIndex[[1]],
        main = "", ylab = ylab, xlab = "",
        ylim = ylims, cex.axis = 1.5, cex.lab = 1.75,
        pch = 21, cex = 1.4, xaxt = "n", bg = pointColor,
        xlim = c(range(age[fIndex[[1]]]) + c(-0.01, 0.07)), ...)
    
    if(smoothIt) {
        lines(age[fIndex[[1]]], fit[fIndex[[1]]], col = lineColor, lwd = 6)
    }
    
    axis(1, at = unique(age[fIndex[[1]]]), 
        labels = 40 + 52 * signif(unique(age[fIndex[[1]]]), 1), cex.axis = 1.5)
    
    if(ageLabel == "bottom") {
        text(x = quantile(age[fIndex[[1]]], 0.33), y = min(ylims), "PCW",
            cex=1.5)
    } else if(ageLabel == "top") {
        text(x = quantile(age[fIndex[[1]]], 0.33), y = max(ylims), "PCW",
            cex=1.5)
    } 
    
    # infant + child
    par(mar = c(4, 0.25,3,0.25))
    for(j in 2:3) {
        plot(y ~ age, subset=fIndex[[j]],
            main = "", ylab = "" , xlab = "", yaxt = "n", cex = 1.4,
            xlim = range(age[fIndex[[j]]]) + c(-0.03, 0.03),
            ylim = ylims, cex.axis = 1.5, pch = 21, bg = pointColor)

        if(ageBreaks[2] == 0 & smoothIt) {
            lines(age[fIndex[[j]]], fit[fIndex[[j]]], col = lineColor, lwd = 6)
        }
        if(ageBreaks[2] < 0 & smoothIt) {
            lines(age[fIndex[[j]]][-1], fit[fIndex[[j]]][-1], col= lineColor,
                lwd = 6)
        }
    }
    
    # adults
    par(mar = c(4, 0.25,3,1))
    plot(y ~ age, subset=fIndex[[4]],
        main = "", ylab = "", xlab = "", yaxt = "n", cex = 1.4,
        xlim = range(age[fIndex[[4]]]) + c(-0.01, 0.01),
        ylim = ylims, cex.axis = 1.5, pch = 21, bg = pointColor)

    if(smoothIt) {
        lines(age[fIndex[[4]]], fit[fIndex[[4]]], col = lineColor, lwd = 6)
    }

    mtext(mainText, outer = TRUE, line = -2.5, cex = 1.35)    
    mtext("Age", side = 1, outer = TRUE, line = -1.5, cex = 1.35)
}