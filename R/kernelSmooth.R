#takes a sparse matrix
smoothHeatmap = function(sm, coords, label=NULL) {

    scale.factor <- 2*log5(20)*log5(100)/log5(ncol(sm)/2) - 2*log5(4)
    message("scale.factor: ", scale.factor)

    # settings no longer options, can experiment later
    # for 1000 sequences this means more bins than sequences !?
    width <- coords[2] - coords[1]
    nbin <- c(round(width*scale.factor), nrow(sm))
    message("nbin: ", nbin[1], ", ", nbin[2])
    bw <- c(3/scale.factor,3)
    message("bw: ", bw[1], ", ", bw[2])

    message("\nCalculating density...")
    df <- summary(sm)
    map <- bkde2D(cbind(df$j, df$i), bandwidth=bw, gridsize=nbin)
    map$fhat <- sum(sm)*map$fhat

    if(is.null(coords)) coords=c(0, ncol(sm))
    new("Heatmap",
        xm=map$x1,
        ym=map$x2,
        value=map$fhat,
        nseq=nrow(sm),
        coords=as.integer(coords),
        max_value=max(map$fhat),
        label=label)
}

log5 <- function(x) {
    log10(x)/log10(5)
}

# heatmap = smoothHeatmap(sm, coords=c(-500, 500), label="TA", transform=function(x) x^(1/3))

