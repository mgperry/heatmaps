#takes a sparse matrix
setGeneric("smooth", function(heatmap, ...) {
    StandardGeneric("smooth")
})

setMethod("smooth", signature(heatmap="Heatmap"),
    function(heatmap, nbin=NULL, vsmooth=1, hsmooth=1, bw=NULL) {
    message("\nCalculating density...")

    if (!is.null(nbin)) {
        if (!length(nbin) == 2) stop("nbin must have length 2")
        if (!all(nbin %% 1 == 0)) stop("nbin must have integer values")
        if (!(vsmooth == 1 && hsmooth == 1)) warning("nbin is set; overriding v/hsmooth")
    } else {
        nbin <- ceiling(c(heatmap@nseq/vsmooth, width(heatmap)/hsmooth))
    }

    if (!all(heatmap@matrix %in% c(0,1))) warning("smooth expects a binary matrix")
    if (all(heatmap@matrix == 0)) warning("no patterns detected")

    message("nbin: ", paste(nbin, collapse=" "))

    if (is.null(bw)) {
        bw <- c(3,3)
    }

    sm <- as(heatmap@matrix, "sparseMatrix")

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

