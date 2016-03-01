plotHeatmapMeta = function(hm_list, binsize, colors, addReferenceLine=FALSE) {
    if (!length(unique(lapply(hm_list, function(x) x@coords))) == 1)
        stop("heatmaps must have the same coordinates")

    if (!length(unique(lapply(hm_list, function(x) x@coords))) == 1)
        stop("heatmaps must have the same number of sequences")

    coords = hm_list[[1]]@coords
    breaks = seq(0, width(hm_list[[1]]), by=binsize)
    bin_sums = lapply(hm_list, bin_heatmap, breaks=breaks)
    occurrence = lapply(bin_sums, function(x) x/(hm_list[[1]]@nseq*binsize))
    max_value = max(vapply(occurrence, max, numeric(1)))
    plot(0, 0, xlim=coords, ylim=c(0, max_value), axes=FALSE, type="n")
    axis(1)
    axis(2)
    for (i in seq_along(occurrence)) {
        x_coord = breaks[1:(length(breaks)-1)] + binsize/2 + coords[1]
        lines(x_coord, occurrence[[i]], col = colors[i], type='l', lwd=2)
    }
}

bin_heatmap = function(hm, breaks) {
    partition = data.table(pos=hm@xm, value=colSums(hm@matrix), bin=cut(hm@xm, breaks))
    partition[, list(sum=sum(value)), by=bin][,sum]
}

plotHeatmapMetaSmooth = function(hm_list, colors, span=0.1, addReferenceLine=FALSE) {
    if (!length(unique(lapply(hm_list, function(x) x@coords))) == 1)
        stop("heatmaps must have the same coordinates")

    if (!length(unique(lapply(hm_list, function(x) x@coords))) == 1)
        stop("heatmaps must have the same number of sequences")

    coords = hm_list[[1]]@coords
    pred = list()
    for (i in seq_along(hm_list)) {
        hm = hm_list[[i]]
        col_sums = colSums(hm@matrix)
        lo = loess(col_sums ~ hm@xm, span=span)
        pred[[i]] = predict(lo)
    }
    max_value = max(vapply(pred, max, numeric(1)))
    plot(0, 0, xlim=coords, ylim=c(0, max_value), axes=FALSE, type="n", xlab=NULL, ylab=NULL)
    axis(1)
    axis(2)
    for (i in seq_along(occurrence)) {
        lines(hm_list[[i]]@xm + coords[1], pred[[i]], col = colors[i], type='l', lwd=2)
    }
}

