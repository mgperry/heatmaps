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

