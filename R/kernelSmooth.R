#takes a sparse matrix
smoothHeatmap = function(sm, coords=NULL, label=NULL, transform=NULL) {

    scale.factor <- 2*log5(20)*log5(100)/log5(ncol(sm)/2) - 2*log5(4)
    message("scale.factor: ", scale.factor)

    # settings no longer options, can experiment later
    # for 1000 sequences this means more bins than sequences !?
    nbin <- c(round(flank*scale.factor), nrow(sm))
    message("nbin: ", nbin[1], ", ", nbin[2])
    bw = c(3/scale.factor,3)
    message("bw: ", bw[1], ", ", bw[2])

    message("\nCalculating density...")
    df = summary(sm)
    map <- bkde2D(cbind(df$i, df$j), bandwidth=bw, gridsize=nbin)

    if(is.null(coords)) coords=c(0, ncol(sm))
    new("Heatmap",
        xm=map$x1,
        ym=map$x2,
        value=map$fhat,
        nseq=nrow(sm),
        coords=as.integer(coords),
        max_value=max(map$fhat)*sum(sm),
        label=label,
        transform=transform)
}

heatmap = smoothHeatmap(sm, coords=c(-500, 500), label="TA", transform=function(x) x^(1/3))

log5 <- function(x) {
    log10(x)/log10(5)
}

