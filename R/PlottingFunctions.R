.smoothScatter <- function(x, y=NULL, map, nbin=128, bandwidth,
colramp=colorRampPalette(c("white", blues9)), nrpoints=100, pch=".", cex=1,
col="black", transformation=function(x) x^.25, max.value=NULL, postPlotHook=box,
xlab=NULL, ylab=NULL, xlim, ylim, xaxs=par("xaxs"), yaxs=par("yaxs"), ...){
    
    if (!is.numeric(nrpoints) | (nrpoints<0) | (length(nrpoints)!=1))
    stop("'nrpoints' should be numeric scalar with value >= 0.")
    
    ## similar as in plot.default
    xlabel <- if (!missing(x)) deparse(substitute(x))
    ylabel <- if (!missing(y)) deparse(substitute(y))
    xy <- xy.coords(x, y, xlabel, ylabel)
    xlab <- if (is.null(xlab)) xy$xlab else xlab
    ylab <- if (is.null(ylab)) xy$ylab else ylab
    
    ## eliminate non-finite (incl. NA) values
    x <- cbind(xy$x, xy$y)[ is.finite(xy$x) & is.finite(xy$y), , drop = FALSE]

    ## xlim and ylim
    if(!missing(xlim)) {
        stopifnot(is.numeric(xlim), length(xlim)==2, is.finite(xlim))
        x <- x[ min(xlim) <= x[,1] & x[,1] <= max(xlim), ]
    } else {
        xlim <- range(x[,1])
    }
    if(!missing(ylim)) {
        stopifnot(is.numeric(ylim), length(ylim)==2, is.finite(ylim))
        x <- x[ min(ylim) <= x[,2] & x[,2] <= max(ylim), ]
    } else {
        ylim <- range(x[,2])
    }
    
    ## create density map [ code in --> ../../grDevices/R/smooth2d.R ]:
    xm <- map$x1
    ym <- map$x2
    dens <- map$fhat
    dens[] <- transformation(dens)

    breaks <- seq(0, transformation(max.value), length.out=257)
    
    ## plot color image
    image(xm, ym, z=dens, col=colramp(256), breaks=breaks, xlab=xlab, ylab=ylab,
        xlim=xlim, ylim=ylim, xaxs=xaxs, yaxs=yaxs, ...)
    if(!is.null(postPlotHook))
    postPlotHook()
    
    ## plot selection of dots
    if (nrpoints > 0) {
        nrpoints <- min(nrow(x), ceiling(nrpoints))
        ## we assume that map$x1 and map$x2 go linearly from
        ## their first to their last value in nbin steps
        stopifnot((nx <- length(xm)) == nrow(dens),
            (ny <- length(ym)) == ncol(dens))
        ixm <- 1L + as.integer((nx-1)*(x[,1]-xm[1])/(xm[nx]-xm[1]))
        iym <- 1L + as.integer((ny-1)*(x[,2]-ym[1])/(ym[ny]-ym[1]))
        sel <- order(dens[cbind(ixm, iym)])[seq_len(nrpoints)]
        points(x[sel,], pch=pch, cex=cex, col=col)
    }
}

.myColorPalette <- function(colorName){
    colors.df <- data.frame(midcol=c("palegreen","aquamarine","lightblue",
        "lavender","mistyrose","peachpuff","lightgoldenrod2","wheat2",
        "lightgray"), highcol=c("darkgreen","darkslategrey","blue4","purple4",
        "deeppink4","red4","darkorange4","salmon4","black"),
        stringsAsFactors = FALSE)
    rownames(colors.df) <- c("green","cyan","blue","purple","pink","red",
        "orange","brown","gray")
    c("white", colors.df[colorName,"midcol"], colors.df[colorName,"highcol"])
}

log5 <- function(x) {
    log10(x)/log10(5)
}

.pattern.smoothscatter <- function(melted, orig, patterns, flankUp=NULL,
flankDown=NULL, bw=NULL, nbin=NULL, color='blue', transf=NULL, xTicks=NULL,
xTicksAt=NULL, xLabel="", yTicks=NULL, yTicksAt=NULL, yLabel="", cex.axis=8,
plot.scale=TRUE, scale.length=NULL, scale.width=10, add.label=TRUE, cex.label=8,
label.col='black', addReferenceLine=TRUE, plotColorLegend=TRUE, out,
plot.width=2000, plot.height=2000, useMulticore=FALSE, nrCores=NULL){
    
    
    mycols <- .myColorPalette(color)

    flank <- width(orig)[1]
    nr.seq <- length(orig)
    
    scale.factor <- 2*log5(20)*log5(100)/log5(flank/2) - 2*log5(4)

    if(length(nbin) == 0){
        nbin <- c(round(flank*scale.factor), length(orig))
    }
    if(length(bw) == 0){
        bw = c(3/scale.factor,3)
    }
    
    sums <- vector()
    for(di in patterns) {
        melted[[di]]$sequence <- length(orig) + 1 - melted[[di]]$sequence
        sums <- c(sums, sum(melted[[di]][,3]))
    }
    names(sums) <- patterns
    
    message("\nCalculating density...")

    if(useMulticore){
        a0 <- mclapply(as.list(patterns), function(di){
            bkde2D(melted[[di]][,c(2,1)], bandwidth=bw, gridsize=nbin);
        }, mc.cores = nrCores)
        names(a0) <- patterns
    }else{
        a0 <- list()
        for(di in patterns){
            message("->", di)
            a0[[di]] <- bkde2D(melted[[di]][,c(2,1)],bandwidth=bw,gridsize=nbin)
        }
    }
    
    d <- vector()
    a <- list()
    for(di in patterns){
        a[[di]] <- a0[[di]]$fhat
        d <- c(d, max(a[[di]]))
    }
    max.value <- max(d*sums)/sums
    names(max.value) <- patterns
    
    if(length(transf) == 0){
        transf <- function(x) {x^(1/3)}
        untransf <- function(x) {x^3}
    }
    
    if(length(xTicks) == 0){
        xTicks <- c(-1*flankUp, -1*flankUp/2, 0, flankDown/2, flankDown)
        xTicksAt <- cumsum(c(0.5, flankUp/2, flankUp/2, flankDown/2-1,
            flankDown/2))
    }
    
    cols.to.draw <- c(0:(flank-1))
    rows.to.draw <- 1:length(orig)
    
    message("\nPlotting...")
    
    if(useMulticore){
        
        mclapply(as.list(patterns), function(di){
            outfile <- paste(out,di, "png", sep=".")
            png(filename=outfile, width = plot.width, height = plot.height)
            par(mar = c(12, 8.5, 2, 8.5))
            dinuc.subset <- melted[[di]]
            
            .smoothScatter(x=dinuc.subset$position-cols.to.draw[1],
            y=max(rows.to.draw)+1-dinuc.subset$sequence, map=a0[[di]],
            axes=FALSE, xaxs="i", yaxs="i", xlab=xLabel, nrpoints=0,
            ylab=yLabel, xlim=c(0.5,flank -0.5), main='', bandwidth=bw,
            nbin=nbin, colramp=colorRampPalette(mycols), transformation=transf,
            max.value=max.value[di], pch=20, cex=0.8, col="gray",
            cex.main=1.5, cex.lab = 1.2)
            
            if(length(xTicks) > 0){
                axis(1, at=xTicksAt, labels=xTicks, cex.axis=cex.axis, padj=1,
                lwd=6, tcl=-1)
            }
            if(length(yTicks) > 0){
                axis(2, at=yTicksAt, labels=yTicks, cex.axis=cex.axis, las=1,
                lwd=6, tcl=-1)
            }
            box(lwd = 6)
            if(plot.scale){
                if(length(scale.length) == 0){
                    scale.length = flank/5
                }
                lines(c(0.03*flank/2,0.03*flank/2 + scale.length),
                    c(0.03*max(dinuc.subset$sequence), 0.03*
                    max(dinuc.subset$sequence)), lwd=scale.width, col=label.col)
                text(x=0.03*flank/2 + scale.length/2,
                    y=0.06*max(dinuc.subset$sequence),
                    labels=paste(round(scale.length), 'bp', sep=''),
                    cex=cex.label, adj=c(0.5,0), col=label.col, font=2)
            }
            if(add.label){
                text(x=0.02*flank/2, y=0.98*max(dinuc.subset$sequence),
                labels=di, cex=cex.label, adj=c(0,1), col=label.col, font=2)
            }
            if(addReferenceLine){
                abline(v=flankUp+0.5, lty="dashed", lwd=6)
            }
            dev.off()
        })
        
    }else{
        
        for (di in patterns){
            message("->", di)
            outfile <- paste(out,di,"png",sep=".")
            png(filename=outfile, width=plot.width, height=plot.height)
            par(mar=c(12, 8.5, 2, 8.5))
            dinuc.subset=melted[[di]]
            
            .smoothScatter(x=dinuc.subset$position-cols.to.draw[1],
            y=max(rows.to.draw)+1-dinuc.subset$sequence, map=a0[[di]],
            axes=FALSE, xaxs="i", yaxs="i", xlab=xLabel,nrpoints=0,
            ylab=yLabel, xlim=c(0.5,flank-0.5), main='', bandwidth=bw,
            nbin=nbin, colramp=colorRampPalette(mycols), transformation=transf,
            max.value=max.value[di], pch=20, cex=0.8,col="gray",
            cex.main=1.5, cex.lab=1.2)
            
            if(length(xTicks) > 0){
                axis(1, at=xTicksAt, labels=xTicks, cex.axis=cex.axis, padj=1,
                lwd=6, tcl=-1)
            }
            if(length(yTicks) > 0){
                axis(2, at=yTicksAt, labels=yTicks, cex.axis=cex.axis, las=1,
                lwd=6, tcl=-1)
            }
            box(lwd = 6)
            if(plot.scale){
                if(length(scale.length) == 0){
                    scale.length = flank/5
                }
                lines(c(0.03*flank/2,0.03*flank/2 + scale.length),
                c(0.03*max(dinuc.subset$sequence), 0.03*
                max(dinuc.subset$sequence)), lwd=scale.width, col=label.col)
                text(x=0.03*flank/2 + scale.length/2,
                y=0.06*max(dinuc.subset$sequence),
                labels=paste(round(scale.length), 'bp', sep=''),
                cex=cex.label, adj=c(0.5,0), col=label.col, font=2)
            }
            if(add.label){
                text(x=0.02*flank/2, y=0.98*max(dinuc.subset$sequence),
                labels=di, cex=cex.label, adj=c(0,1), col=label.col, font=2)
            }
            if(addReferenceLine){
                abline(v=flankUp+0.5, lty="dashed", lwd=6)
            }
            dev.off()
        }
    }
    
    if(plotColorLegend == TRUE){
        png(filename=paste(out,"ColorLegend","png",sep="."),
        width=0.15*plot.height, height=plot.height)
        f <- colorRampPalette(mycols)
        nr.labels <- 4
        leg <- rep('', 256)
        leg[seq(0.05*256, 0.95*256, length.out=nr.labels)] <-
            formatC(untransf(seq(0.05*transf(max(d*sums)), 0.95*
            transf(max(d*sums)), length.out=nr.labels)), format='f', digits=2)
        par(mar = c(12, 14, 2, 0.5))
        plot(1, 1, type='n', bty='n', xaxt='n', yaxt='n', xlim=c(0,1),
        ylim=c(0,7), xaxs="i", yaxs="i", xlab='', ylab='')
        box(lwd = 6)
        color.legend(0, 0, 1, 7, legend=leg, rect.col=f(256), align="lt",
        gradient='y', cex=cex.axis)
        dev.off()
    }
    
    return(a)
}


#################################
# Function for plotting heatmap


.plot.motif.heatmap <- function(motifScanningScores, flankUp=NULL,
flankDown=NULL, cols, breaks, xTicks=NULL, xTicksAt=NULL, xLabel="",
yTicks=NULL, yTicksAt=NULL, yLabel="", cexAxis=8, plotScale=TRUE,
scaleLength=NULL, scaleWidth=15, addReferenceLine=TRUE){
    
    par(mar = c(12, 8.5, 2, 8.5))
    
    if(length(xTicks) == 0){
        xTicks <- c(-1*flankUp, -1*flankUp/2, 0, flankDown/2, flankDown)
        xTicksAt <- cumsum(c(0.5, flankUp/2, flankUp/2, flankDown/2-1,
            flankDown/2))
    }
    
    image(t(motifScanningScores[nrow(motifScanningScores):1,]), col=cols,
    breaks=breaks, xaxt="n", yaxt="n")
    
    if(length(xTicks) > 0){
        axis(1, at=xTicksAt/(flankUp+flankDown), labels=xTicks,
        cex.axis=cexAxis, padj=1, lwd=6, tcl=-1)
    }
    if(length(yTicks) > 0){
        axis(2, at=yTicksAt/nrow(motifScanningScores), labels=yTicks,
        cex.axis=cexAxis, las=1, lwd=6, tcl=-1)
    }
    box(lwd = 6)
    if(plotScale){
        if(length(scaleLength) == 0){
            scaleLength <- (flankUp + flankDown)/5
        }
        lines(c(0.015, 0.015+scaleLength/(flankUp + flankDown)),
        c(0.03,0.03), lwd=scaleWidth, col="gray90")
        text(x=0.03/2+scaleLength/(2*(flankUp + flankDown)), y=0.06,
        labels=paste(round(scaleLength), 'bp', sep=''), cex=cexAxis,
        adj=c(0.5,0), col="gray90", font=2)
    }
    if(addReferenceLine){
        abline(v=(flankUp+0.5)/(flankUp+flankDown), lty="dashed",
        col="gray90", lwd=6)
    }

}



#######################################
# Function of plotting average signal


.plot.windowed.average <- function(occurence.melted.list, nr.seq,
pattern.widths, flankUp = NULL, flankDown = NULL, smoothingWindow = 3,
color = rainbow(length(occurence.melted.list)), xLabel =
"Distance to reference point (bp)", yLabel = "Relative frequency", cexAxis = 2,
addReferenceLine = TRUE, plotLegend = TRUE, cexLegend = 2, add = FALSE, ...){
    
    message("\nPlotting average signal...\n")
    patterns.avg.signal<-lapply(c(1:length(occurence.melted.list)),
    function(i) {
        x <- occurence.melted.list[[i]]
        a.s <- table(x$position)/nr.seq
        avg.signal <- rep(0, times =
        flankUp+flankDown-pattern.widths[i]+1)
        names(avg.signal) <- c(1:length(avg.signal))
        avg.signal[names(a.s)] <- a.s
        return(avg.signal)
    })
    
    starts <- lapply(patterns.avg.signal, function(x) {
        seq(1, length(x), by = smoothingWindow)
    })
    ends <- lapply(starts, function(x) {x + smoothingWindow - 1})
    ends <- lapply(1:length(ends), function(x) {
        ends[[x]][length(ends[[x]])] <- min(ends[[x]][length(ends[[x]])],
        length(patterns.avg.signal[[x]])); return(ends[[x]])
    })
    
    patterns.avg.signal.windowed <- lapply(c(1:length(patterns.avg.signal)),
    function(x){
        a.s <- sapply(seq(1:length(starts[[x]])), function(y) {
            mean(patterns.avg.signal[[x]][starts[[x]][y]:ends[[x]][y]])})
        x.coor <- (starts[[x]]+ends[[x]])/2 - flankUp - 1
        return(list(x.coor, a.s))
    })
    
    if(!(add)){
        plot(1, 1, xlim = c(-flankUp, flankDown), ylim = c(0,
        1.05*max(unlist(lapply(patterns.avg.signal.windowed, function(x) {
            max(x[[2]])})))), type = "n", xlab = xLabel, ylab = yLabel,
        cex.axis = cexAxis, cex.lab = cexAxis, ...)
    }
    
    a.s <- lapply(c(1:length(patterns.avg.signal.windowed)), function(i){
        lines(x = patterns.avg.signal.windowed[[i]][[1]],
        y = patterns.avg.signal.windowed[[i]][[2]], type = "l",
        col = color[i], ...)
    })
    
    if(plotLegend){
        legend("top", legend = names(pattern.widths), bty = "n", horiz = TRUE,
        lwd = 1, col = color, cex = cexLegend)
    }
    if(addReferenceLine){
        abline(v = 0, lty = "dashed")
    }


}

