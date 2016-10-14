#' Plot a Meta-region plot from heatmaps
#'
#' @param hm_list A list of heatmaps
#' @param binsize Integer, size of bins to use in plot
#' @param colors Color to use for each heatmap
#' @param addReferenceLine Logical, add reference line at zero or not
#'
#' This function creates a meta-region plot from 1 or more heatmaps with the same
#' coordinates. A meta-region plot graphs the sum of the signal at each position in
#' each heatmap rather than visualising the signal in two dimensions. Often binning
#' is required to smooth noisy signal.
#'
#' @return invisible(0)
#'
#' @export
#' @importFrom graphics plot mtext legend axis lines
#' @examples
#' data(HeatmapExamples)
#' plotHeatmapMeta(hm, color="steelblue")
plotHeatmapMeta = function(hm_list, binsize=1, colors=gg_col(length(hm_list)), addReferenceLine=FALSE) {
    if (class(hm_list) == "Heatmap") hm_list = list(hm_list) # allow single heatmap argument

    if (!length(unique(lapply(hm_list, function(x) x@coords))) == 1)
        stop("heatmaps must have the same coordinates")

    co_ord = coords(hm_list[[1]])

    if (!length(unique(lapply(hm_list, function(x) xm(x)))) == 1)
        stop("heatmaps must have the same xm values")

    n_seq = unique(sapply(hm_list, function(x) x@nseq))

    if (!length(n_seq) == 1)
        stop("heatmaps must have the same number of sequences")

    coords = hm_list[[1]]@coords
    if (binsize != 1) {
        if (!all(xm(hm_list[[1]]) == 1:width(hm_list[[1]]))) {
            stop("cannot set binsize for heatmaps which are already binned/smoothed")
        }
        breaks = seq(0, width(hm_list[[1]]), by=binsize)
        bin_sums = lapply(hm_list, bin_heatmap, breaks=breaks)
        x_coord = breaks[1:(length(breaks)-1)] + binsize/2 + coords[1]
    } else {
        breaks = xm(hm_list[[1]])
        x_coord = breaks + binsize/2 + coords[1]
        bin_sums = lapply(hm_list, function(x) colSums(image(x)))
    }
    scale_factor = n_seq*(width(hm_list[[1]])/length(x_coord))
    occurrence = lapply(bin_sums, function(x) x/scale_factor)
    max_value = max(vapply(occurrence, max, numeric(1)))
    plot(0, 0, xlim=coords, ylim=c(0, max_value), axes=FALSE, type="n", xlab="", ylab="")
    axis(1)
    axis(2)
    mtext("Relative Position", side = 1, line = 3, cex = 1.5, font = 2)
    mtext("Frequency", side = 2, line = 2.5, cex = 1.5, font = 2)
    for (i in seq_along(occurrence)) {
        lines(x_coord, occurrence[[i]], col = colors[i], type='l', lwd=2)
    }
    labels = vapply(hm_list, function(x) x@label, character(1))
    legend('topright', labels, bty="n", lty=1, lwd=2, col=colors)
}

#' @importFrom stats aggregate
bin_heatmap = function(hm, breaks) {
    partition = data.frame(pos=xm(hm), value=colSums(image(hm)), bin=cut(xm(hm), breaks))
    aggregate(partition$value, sum, by=list(bin=partition$bin))$x
}

#' @importFrom grDevices hcl
gg_col = function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
}


