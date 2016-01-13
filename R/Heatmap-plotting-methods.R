setMethod("plot", signature="Heatmap",
    function(heatmap, options=NULL, ...) {

    if (is.null(options)) {
        options = heatmapOptions(...)
    }

    colramp = .myColorPalette(options$color) # could check for function as input

    breaks <- seq(0, heatmap@transform(heatmap@max_value), length.out=257) # could be option

    xm <- heatmap@xm
    ym <- heatmap@ym
    val <- options$transform(heatmap@value)

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
        lenged.width=0.2,
        addReferenceLine=TRUE,
        transform=function(x) x,
        untransform=function(x) x)
    def[names(usr)] = usr
    return(def)
}

plotHeatmapList = function(heatmap_list, groups=NULL, group_colours=NULL, group_transform=NULL, options) {
    n_plots = length(heatmap_list)

    # colors for groups?

    if (is.null(groups) || length(unique(groups)) == n_plots) {
        groups = 1:n_plots
    } else if (length(groups) == 1) {
        if (groups != 1) stop("groups incorrectly formatted, refer to ?plotHeatmapList")
        groups = rep(1, n_plots)
    } else if (length(groups) != n_plots) {
        stop("groups must have 1 value per plot")
    } else {
        groups = as.numeric(factor(groups, levels=unique(groups)))
    }

    if (groups != 1:n_plots) {
        for (i = 1:max(groups)) {
            group = which(groups == i)
            max_d = sum(vapply(heatmaps_list[group], function(x) max(x@value), integer(1)))
            for (index in group) {
                heatmaps_list[[index]]@max_value = max_d
            }
        }
    }

    group_list = split(1:n_plots, groups)

    if (options$legend = TRUE) {
        widths = numeric(0)
        for (g in group_list) {
            widths = c(widths, 0.2, rep(1, length(g)))
        }
        mat = 1:length(widths)
        widths = widths/sum(widths)
        layout(mat, widths)
    } else {
        layout(1:7)
    }

    for (g in group_list) {
        if(options$legend = TRUE) {
            plot_legend(heatmap_list[g[1]], options)
        }
        for (i in g) {
            plot(heatmap_list[i], options)
        }
    }
}


# additional functions

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
    colnames = c("green", "cyan", "blue", "purple", "pink", "red", "orange", "brown", "gray")
        if(!(colorName %in% colnames)){
            msg = paste("Colour", colorName, "not supported, must be one of", paste(colnames, collapse=', '))
            stop(msg)
        }
    colors.df <- data.frame(
        midcol=c("palegreen", "aquamarine", "lightblue", "lavender",
                 "mistyrose", "peachpuff", "lightgoldenrod2", "wheat2", "lightgray"),
        highcol=c("darkgreen", "darkslategrey", "blue4", "purple4",
                  "deeppink4", "red4", "darkorange4", "salmon4", "black"),
        stringsAsFactors = FALSE)
    rownames(colors.df) <- colnames
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

