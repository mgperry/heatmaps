#takes a sparse matrix
setGeneric("smooth", function(heatmap, ...) {
    StandardGeneric("smooth")
})

setMethod("smooth", signature(heatmap="Heatmap"),
    function(heatmap, nbin=NULL, bw=NULL) {

    scale.factor <- 2*log5(20)*log5(100)/log5(ncol(heatmap@matrix)/2) - 2*log5(4)
    message("scale.factor: ", scale.factor)

    # settings no longer options, can experiment later
    # for 1000 sequences this means more bins than sequences !?
    width <- heatmap@coords[2] - heatmap@coords[1]
    if (is.null(nbin)) {
        nbin <- c(nrow(heatmap@matrix), round(width*scale.factor))
        message("nbin: ", nbin[1], ", ", nbin[2])
    } else {
        if (!length(nbin) == 2) stop("nbin must have length 2")
        if (!all(nbin %% 1 == 0)) stop("nbin must have integer values")
    }
    if (is.null(bw)) {
        bw <- c(3,3/scale.factor)
        message("bw: ", bw[1], ", ", bw[2])
    }

    message("\nCalculating density...")

    sm <- as(heatmap@matrix, "sparseMatrix")
    total_value = sum(sm)

    df <- summary(sm)
    map <- bkde2D(cbind(df$i, df$j), bandwidth=bw, gridsize=nbin, range.x=list(c(1, nrow(heatmap@matrix)), c(1, width)))
    map$fhat <- sum(heatmap@matrix)*map$fhat
    heatmap@xm = map$x2
    heatmap@ym = map$x1
    heatmap@matrix = map$fhat
    heatmap@max_value = max(map$fhat)
    return(heatmap)
})

log5 <- function(x) {
    log10(x)/log10(5)
}

