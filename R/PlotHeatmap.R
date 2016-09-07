setGeneric(name="plotHeatmap",
           def=function(heatmap, options=NULL, ...) {
    standardGeneric("plotHeatmap")
})

setMethod("plotHeatmap", signature="Heatmap",
    function(heatmap, options=NULL, ...) {

    if (is.null(options)) {
        options = heatmapOptions(...)
    }

    if (length(options$color) == 1) {
        col_palette = colorRampPalette(default_color(options$color))
    } else {
        col_palette = colorRampPalette(options$color)
    }

    if (!is.na(options$transform)) {
        val = options$transform(heatmap@matrix)
        breaks = seq(0, options$transform(heatmap@max_value), length.out=257)
    } else {
        val = heatmap@matrix
        breaks = seq(0, heatmap@max_value, length.out=257)
    }
    val[val < 0] = 0
    val[val > heatmap@max_value] = heatmap@max_value

    image(xm(heatmap), ym(heatmap), z=t(val),
          col=col_palette(256), breaks=breaks,
          xlim=c(0.5,width(heatmap)-0.5), # remove fluffy edges
          xlab="", ylab="",
          axes=FALSE)

    if (options$box.width > 0) box(lwd = options$box.width)

    if(options$x.ticks == TRUE) {
        x.ticks = make_x_ticks(heatmap@coords)
        axis(1, at=x.ticks, labels=names(x.ticks), cex.axis=options$cex.axis, lwd=options$box.width, tcl=options$tcl, padj=options$padj)
    }

    if(options$scale){
        scale.length = width(heatmap)/5
        x.start.in = grconvertX(0, from="user", to="inches") + options$label.xpos
        x.start = grconvertX(x.start.in, from="inches", to="user")
        y.start.in = grconvertY(0, from="user", to="inches") + options$label.ypos
        y.start = grconvertY(y.start.in, from="inches", to="user")
        lines(c(x.start, x.start + scale.length),
              c(y.start, y.start),
              lwd=options$scale.lwd, col=options$label.col)
        text(x=x.start, y=y.start*2,
             labels=paste(round(scale.length), 'bp', sep=''),
             cex=options$cex.scale, adj=c(0,0), col=options$label.col, font=2)
    }

    if(options$label){
        x.start.in = grconvertX(0, from="user", to="inches") + options$label.xpos
        x.start = grconvertX(x.start.in, from="inches", to="user")
        y.start.in = grconvertY(0, from="user", to="inches") + options$label.ypos
        y.start = grconvertY(y.start.in, from="inches", to="user")
        text(x=x.start, y=length(heatmap)-y.start,
        labels=heatmap@label, cex=options$cex.label, adj=c(0,1), col=options$label.col, font=2)
    }

    if(options$refline){
        abline(v=(-heatmap@coords[1])+0.5, lty="dashed", lwd=options$refline.width)
    }
})

heatmapOptions = function(...) {
    opts = list(...)
    default = list(
        color='Blues',
        box.width=1.25,
        x.ticks=TRUE,
        tcl=-0.25,
        padj=-1.8,
        scale=FALSE,
        scale.lwd=3,
        label.xpos=0.1,
        label.ypos=0.1,
        cex.axis=0.8,
        label=TRUE,
        cex.label=1.2,
        label.col='black',
        legend=FALSE,
        legend.width=0.2,
        legend.ticks=5,
        legend.pos='l',
        cex.legend=1,
        refline=TRUE,
        refline.width=0.5,
        transform=NA,
        plot.mai=list(c(0.3, 0.15, 0.05, 0.15)),
        legend.mai=list(c(0.3, 0.3, 0.05, 0.01))
    )
    default[names(opts)] = opts
    return(default)
}

make_x_ticks = function(coord) {
    xTicks = c(1*coord[1], "", 0, "", coord[2])
    xTicksAt = cumsum(c(0.5, -coord[1]/2, -coord[1]/2, coord[2]/2-1, coord[2]/2))
    names(xTicksAt) = xTicks
    xTicksAt
}

default_color = function(col) {
    if (col == "rainbow") {
        palette = rev(rainbow(9, start=0, end=4/6))
    } else {
       palette = tryCatch(brewer.pal(9, col), error=color_error)
    }
    palette
}

color_error = function(e) "Incorrect color specification: refer to ?default_color"


