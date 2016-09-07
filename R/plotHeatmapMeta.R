plotHeatmapMeta = function(hm_list, binsize=1, colors=gg_col(length(hm_list)), addReferenceLine=FALSE) {
    if (!length(unique(lapply(hm_list, function(x) x@coords))) == 1)
        stop("heatmaps must have the same coordinates")

    if (!length(unique(lapply(hm_list, function(x) x@coords))) == 1)
        stop("heatmaps must have the same number of sequences")

    coords = hm_list[[1]]@coords
    breaks = seq(0, width(hm_list[[1]]), by=binsize)
    bin_sums = lapply(hm_list, bin_heatmap, breaks=breaks)
    occurrence = lapply(bin_sums, function(x) x/(hm_list[[1]]@nseq*binsize))
    max_value = max(vapply(occurrence, max, numeric(1)))
    plot(0, 0, xlim=coords, ylim=c(0, max_value), axes=FALSE, type="n", xlab="", ylab="")
    axis(1)
    axis(2)
    mtext("Relative Position", side = 1, line = 3, cex = 1.5, font = 2)
    mtext("Frequency", side = 2, line = 2.5, cex = 1.5, font = 2)
    for (i in seq_along(occurrence)) {
        x_coord = breaks[1:(length(breaks)-1)] + binsize/2 + coords[1]
        lines(x_coord, occurrence[[i]], col = colors[i], type='l', lwd=2)
    }
    labels = vapply(hm_list, function(x) x@label, character(1))
    legend('topright', labels, bty="n", lty=1, lwd=2, col=colors)
}

bin_heatmap = function(hm, breaks) {
    partition = data.frame(pos=xm(hm), value=colSums(hm@matrix), bin=cut(xm(hm), breaks))
    partition[, list(sum=sum(value)), by=bin][,sum]
}

plotHeatmapMetaSmooth = function(hm_list, span=0.1, colors=gg_col(length(hm_list)), addReferenceLine=FALSE) {
    if (!length(unique(lapply(hm_list, function(x) x@coords))) == 1)
        stop("heatmaps must have the same coordinates")

    if (!length(unique(lapply(hm_list, function(x) x@coords))) == 1)
        stop("heatmaps must have the same number of sequences")

    coords = hm_list[[1]]@coords
    pred = list()
    for (i in seq_along(hm_list)) {
        hm = hm_list[[i]]
        col_means = colSums(hm@matrix)/hm@nseq
        lo = loess(col_means ~ xm(hm), span=span)
        pred[[i]] = predict(lo)
    }
    max_value = max(vapply(pred, max, numeric(1)))
    plot(0, 0, xlim=coords, ylim=c(0, max_value), axes=FALSE, type="n", xlab=NULL, ylab=NULL)
    axis(1)
    axis(2)
    for (i in seq_along(occurrence)) {
        lines(xm(hm_list[[i]]) + coords[1], pred[[i]], col = colors[i], type='l', lwd=2)
    }
    mtext("Relative Position", side = 1, line = 3, cex = 1.5, font = 2)
    mtext("Frequency", side = 2, line = 2.5, cex = 1.5, font = 2)
    labels = vapply(hm_list, function(x) x@label, character(1))
    legend('topright', labels, bty="n", lty=1, lwd=2, col=colors)
}

gg_col = function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
}


