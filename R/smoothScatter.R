.smoothScatter <- function(
        map,
        color='blue',
        transf=NULL, 
        max_value=NULL,
        xTicks=xTicks, cex.axis=cex.axis,
        plot.scale=plot.scale, scale.length=scale.length, scale.width=scale.width, flank=flank,
        add.label=add.label, cex.label=cex.label, label.col=label.col, di,
        addReferenceLine=addReferenceLine, flankUp=flankUp, max_dinuc=max_dinuc) {

    if (is.null(transf)) transf = function(x) x
    colramp = .myColorPalette(color)

    if (is.null(max_value)) max_value = max(map$fhat)
    breaks <- seq(0, transf(max_value), length.out=257)

    xm <- map$x1
    ym <- map$x2
    dens <- transf(map$fhat)

    message("plotting in .smoothScatter")
    ## plot color image
    image(xm, ym, z=dens,
          col=colramp(256), breaks=breaks,
          xlim=c(0.5,flank-0.5), # remove fluffy edges
          xaxs="i", yaxs="i", # may be default
          xlab="", ylab="",
          axes=FALSE, # maybe goes here
          # previously passed via dots
          pch=20, cex=0.8,
          main='', cex.main=1.5)
    box(lwd = 6)
    message("plotting labels in .smoothScatter")
    if(length(xTicks) > 0){
        axis(1, at=xTicks, labels=names(xTicks), cex.axis=cex.axis, padj=1,
        lwd=6, tcl=-1)
    }
    if(plot.scale){
        if(length(scale.length) == 0){
            scale.length = flank/5
        }
        lines(c(0.03*flank/2,0.03*flank/2 + scale.length),
                c(0.03*max_dinuc, 0.03*
                max_dinuc), lwd=scale.width, col=label.col)
        text(x=0.03*flank/2 + scale.length/2,
                y=0.06*max_dinuc,
                labels=paste(round(scale.length), 'bp', sep=''),
                cex=cex.label, adj=c(0.5,0), col=label.col, font=2)
    }
    if(add.label){
        text(x=0.02*flank/2, y=0.98*max_dinuc,
        labels=di, cex=cex.label, adj=c(0,1), col=label.col, font=2)
    }
    if(addReferenceLine){
        abline(v=flankUp+0.5, lty="dashed", lwd=6)
    }
}

.pattern.smoothscatter <- function(
    melted, coord, nseq, patterns, out,
    bw=NULL, color='blue', transf=NULL,
    xTicks=NULL, xTicksAt=NULL, yTicks=NULL, yTicksAt=NULL, cex.axis=8,
    plot.scale=TRUE, scale.length=NULL, scale.width=10,
    add.label=TRUE, cex.label=8, label.col='black',
    addReferenceLine=TRUE, plotColorLegend=TRUE) {

    flank <- coord[2] - coord[1]

    scale.factor <- 2*log5(20)*log5(100)/log5(flank/2) - 2*log5(4)
    message("scale.factor: ", scale.factor)

    # settings no longer options, can experiment later
    # for 1000 sequences this means more bins than sequences !?
    nbin <- c(round(flank*scale.factor), nseq)
    message("nbin: ", nbin[1], ", ", nbin[2])
    bw = c(3/scale.factor,3)
    message("bw: ", bw[1], ", ", bw[2])

    sums <- vector()
    for(di in patterns) {
        melted[[di]]$sequence <- nseq + 1 - melted[[di]]$sequence # inverts
        sums <- c(sums, sum(melted[[di]][,3])) # basically equal to length(orig)
    }
    names(sums) <- patterns

    message("\nCalculating density...")

    smoothed_map <- list()
    for(di in patterns){
        message("->", di)
        smoothed_map[[di]] <- bkde2D(melted[[di]][,c(2,1)],bandwidth=bw,gridsize=nbin)
    }

    # max_density <- vector()
    # for(di in patterns){
    #     max_density <- c(max_density, max(smoothed_map[[di]]$fhat))
    # }

    max_density = vapply(patterns, function(di) max(smoothed_map[[di]]$fhat), numeric(1))
    max_value <- max(max_density*sums)/sums
    names(max_value) <- patterns

    message("Max value:")
    message(paste(names(max_value), collapse=" "))
    message(paste(max_value, collapse=" "))

    if(length(transf) == 0){
        transf <- function(x) {x^(1/3)}
        untransf <- function(x) {x^3}
    }

    xTicks = make_x_ticks(coord)
    message("\nPlotting...")

    for (di in patterns){
        message("->", di)
        outfile <- paste(out,di,"png",sep=".")
        png(filename=outfile, width=2000, height=2000)
        par(mar=c(12, 8.5, 2, 8.5))

        max_dinuc = max(melted[[di]]$sequence)
        flankUp=-coord[1]

        .smoothScatter(
            map=smoothed_map[[di]],
            color='blue', 
            transf=transf,
            max_value=max_value[di],
            # scale/label options
            xTicks=xTicks, cex.axis=cex.axis,
            plot.scale=plot.scale, scale.length=scale.length, scale.width=scale.width, flank=flank, 
            add.label=add.label, cex.label=cex.label, label.col=label.col, di,
            addReferenceLine=addReferenceLine, flankUp=flankUp, max_dinuc=max_dinuc) 

        dev.off()
    }

    message("plotting legend")
    max_value_all = transf(max(max_density*sums))
    if(plotColorLegend == TRUE){
        png(filename=paste(out,"ColorLegend","png",sep="."), width=0.15*2000, height=2000)
        plot_legend(max_value_all, untransf, color='blue', cex=8)
        dev.off()
    }

    return(a)
}

plot_legend <- function(max_value, untransf, color='blue', cex=8) {
        f <- .myColorPalette(color)
        nr.labels <- 4
        leg <- rep('', 256)
        leg[seq(0.05*256, 0.95*256, length.out=nr.labels)] <-
            formatC(untransf(seq(0.05*max_value, 0.95*max_value, length.out=nr.labels)), 
                    format='f', digits=2)
        par(mar = c(12, 14, 2, 0.5))
        plot(1, 1, type='n', bty='n', xaxt='n', yaxt='n', xlim=c(0,1),
        ylim=c(0,7), xaxs="i", yaxs="i", xlab='', ylab='')
        box(lwd = 6)
        color.legend(0, 0, 1, 7, legend=leg, rect.col=f(256), align="lt", gradient='y', cex=cex)
}

.myColorPalette <- function(colorName){
    colors.df <- data.frame(
        midcol=c("palegreen", "aquamarine", "lightblue", "lavender", "mistyrose", "peachpuff", "lightgoldenrod2", "wheat2", "lightgray"),
        highcol=c("darkgreen", "darkslategrey", "blue4", "purple4", "deeppink4", "red4", "darkorange4", "salmon4", "black"),
        stringsAsFactors = FALSE)
    rownames(colors.df) <- c("green", "cyan", "blue", "purple", "pink", "red", "orange", "brown", "gray")
    colorRampPalette( c("white", colors.df[colorName,"midcol"], colors.df[colorName,"highcol"]) )
}

log5 <- function(x) {
    log10(x)/log10(5)
}

make_x_ticks = function(coord) {
    xTicks <- c(1*coord[1], 1*coord[1]/2, 0, coord[2]/2, coord[2])
    xTicksAt <- cumsum(c(0.5, -coord[1]/2, -coord[1]/2, coord[2]/2-1,
        coord[2]/2))
    names(xTicksAt) = xTicks
    return(xTicksAt)
}

