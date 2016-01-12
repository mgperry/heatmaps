# should be generic called plotHeatmap or similar
ss_heatmap <- function(
        heatmap,
        options=NULL,
        ...) {

    if (is.null(options)) {
        options = heatmapOptions(...)
    }

    colramp = .myColorPalette(options$color) # could check for function as input

    breaks <- seq(0, heatmap@transform(heatmap@max_value), length.out=257) # could be option

    xm <- heatmap@xm
    ym <- heatmap@ym
    val <- heatmap@transform(heatmap@value)

    message("plotting in ss_heatmap")
    ## plot color image
    image(xm, ym, z=val,
          col=colramp(256), breaks=breaks,
          xlim=c(0.5,width(heatmap)-0.5), # remove fluffy edges
          xaxs="i", yaxs="i", # may be default
          xlab="", ylab="",
          axes=FALSE, # maybe goes here
          # previously passed via dots
          pch=20, cex=0.8,
          main='', cex.main=1.5)

    message("plotting labels in ss_heatmap")

    if (options$box.width > 0) box(lwd = options$box.width)

    if(options$x_ticks == TRUE) {
        x_ticks = make_x_ticks(heatmap@coords)
        axis(1, at=x_ticks, labels=names(x_ticks), cex.axis=options$cex.axis, padj=1, lwd=6, tcl=-1)
    }

    if(options$scale){
        scale.length = width(heatmap)/5
        lines(c(0.03*width(heatmap)/2,0.03*width(heatmap)/2 + scale.length),
                c(0.03*length(heatmap), 0.03*
                length(heatmap)), lwd=options$scale.width, col=options$label.col)
        text(x=0.03*width(heatmap)/2 + scale.length/2,
                y=0.06*length(heatmap),
                labels=paste(round(scale.length), 'bp', sep=''),
                cex=options$cex.label, adj=c(0.5,0), col=options$label.col, font=2)
    }

    if(options$label){
        text(x=0.02*width(heatmap)/2, y=0.98*length(heatmap),
        labels=heatmap@label, cex=options$cex.label, adj=c(0,1), col=options$label.col, font=2)
    }

    if(options$addReferenceLine){
        abline(v=(-heatmap@coords[1])+0.5, lty="dashed", lwd=6)
    }

    # plot legend
}

heatmapOptions = function(...) {
    usr = list(...)
    def = list(
        color='blue',
        box.width=6,
        x_ticks=TRUE,
        cex.axis=8,
        scale=TRUE,
        scale.width=10,
        label=TRUE,
        cex.label=8,
        label.col='black',
        legend=TRUE,
        addReferenceLine=TRUE)
    def[names(usr)] = usr
    return(def)
}

