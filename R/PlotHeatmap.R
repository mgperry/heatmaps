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

    if(options$x_ticks == TRUE) {
        x_ticks = make_x_ticks(heatmap@coords)
        axis(1, at=x_ticks, labels=names(x_ticks), cex.axis=options$cex.axis, padj=1, lwd=6, tcl=-1)
    }

    if(options$scale){
        scale.length = width(heatmap)/5
        x.start.in = grconvertX(0, from="user", to="inches") + 0.8
        x.start = grconvertX(x.start.in, from="inches", to="user")
        y.start.in = grconvertY(0, from="user", to="inches") + 0.8
        y.start = grconvertY(y.start.in, from="inches", to="user")
        lines(c(x.start, x.start + scale.length),
              c(y.start, y.start),
              lwd=options$scale.lwd, col=options$label.col)
        text(x=x.start+scale.length/2, y=y.start*2,
             labels=paste(round(scale.length), 'bp', sep=''),
             cex=options$cex.label, adj=c(0.5,0), col=options$label.col, font=2)
    }

    if(options$label){
        x.start.in = grconvertX(0, from="user", to="inches") + 0.8
        x.start = grconvertX(x.start.in, from="inches", to="user")
        y.start.in = grconvertY(0, from="user", to="inches") + 0.8
        y.start = grconvertY(y.start.in, from="inches", to="user")
        text(x=x.start, y=length(heatmap)-y.start,
        labels=heatmap@label, cex=options$cex.label, adj=c(0,1), col=options$label.col, font=2)
    }

    if(options$addReferenceLine){
        abline(v=(-heatmap@coords[1])+0.5, lty="dashed", lwd=6)
    }
})

heatmapOptions = function(...) {
    opts = list(...)
    def = list(
        color='blue',
        box.width=6,
        x_ticks=TRUE,
        scale=TRUE,
        scale.lwd=15,
        cex.axis=6,
        label=TRUE,
        cex.label=8,
        label.col='black',
        legend=TRUE,
        legend.width=0.2,
        legend.ticks=5,
        legend.pos='l',
        cex.legend=6,
        addReferenceLine=TRUE,
        transform=NA)
    def[names(opts)] = opts
    return(def)
}

default_color <- function(color_name){
    colors = c("green", "cyan", "blue", "purple", "pink", "red", "orange", "brown", "gray")
        if(!(color_name %in% colors)){
            msg = paste("Colour", color_name, "not supported, must be one of", paste(colors, collapse=', '))
            stop(msg)
        }
    colors.df <- data.frame(
        midcol=c("palegreen", "aquamarine", "lightblue", "lavender",
                 "mistyrose", "peachpuff", "lightgoldenrod2", "wheat2", "lightgray"),
        highcol=c("darkgreen", "darkslategrey", "blue4", "purple4",
                  "deeppink4", "red4", "darkorange4", "salmon4", "black"),
        stringsAsFactors = FALSE)
    rownames(colors.df) <- colors
    return(c("white", colors.df[color_name,]))
}

make_x_ticks = function(coord) {
    xTicks <- c(1*coord[1], 1*coord[1]/2, 0, coord[2]/2, coord[2])
    xTicksAt <- cumsum(c(0.5, -coord[1]/2, -coord[1]/2, coord[2]/2-1,
        coord[2]/2))
    names(xTicksAt) = xTicks
    return(xTicksAt)
}

