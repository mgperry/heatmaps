#takes a sparse matrix
setGeneric("smooth", function(heatmap, ...) {
    StandardGeneric("smooth")
})

setMethod("smooth", signature(heatmap="Heatmap"),
    function(heatmap, nbin=NULL, bw=NULL, vsmooth=1, hsmooth=1) {
    message("\nCalculating density...")

    if (!is.null(nbin)) {
        if (!length(nbin) == 2) stop("nbin must have length 2")
        if (!all(nbin %% 1 == 0)) stop("nbin must have integer values")
        if (!(vsmooth == 1 && hsmooth == 1)) warn("nbin is set; overriding v/hsmooth")
    } else {
        nbin <- ceiling(c(heatmap@nseq/vsmooth, width(heatmap)/hsmooth))
    }

    if (is.null(bw)) {
        bw <- c(3,3)
    }

    sm <- as(heatmap@matrix, "sparseMatrix")
    total_value = sum(sm)

    df <- summary(sm)
    map <- bkde2D(cbind(df$i, df$j), bandwidth=bw, gridsize=nbin,
                range.x=list(c(1, heatmap@nseq), c(1, width(heatmap))))

    fhat <- sum(heatmap@matrix)*map$fhat
    heatmap@xm = map$x2
    heatmap@ym = map$x1
    heatmap@matrix = fhat
    heatmap@max_value = max(fhat)
    return(heatmap)
})

