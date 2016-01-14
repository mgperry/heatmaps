setGeneric(name="plotHeatmap",
           def=function(heatmap, options=NULL, ...) {
    standardGeneric("plotHeatmap")
})

setMethod("plotHeatmap", signature="Heatmap",
    function(heatmap, options=NULL, ...) {

    if (is.null(options)) {
        options = heatmapOptions(...)
    }

    colramp = colorRampPalette(.myColorPalette(options$color)) # could check for function as input

    breaks <- seq(0, options$transform(heatmap@max_value), length.out=257) # could be option

    xm <- heatmap@xm
    ym <- heatmap@ym
    val <- options$transform(heatmap@value)

    image(xm, ym, z=val,
          col=colramp(256), breaks=breaks,
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
})

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
        legend.width=0.2,
        legend.ticks=5,
        cex.legend=6,
        addReferenceLine=TRUE,
        transform=function(x) x)
    def[names(usr)] = usr
    return(def)
}

plotHeatmapList = function(heatmap_list, groups=NULL, options=heatmapOptions(), ...) {
    if (class(heatmap_list) == "Heatmap") heatmap_list = list(heatmap_list) # allow single heatmap argument
    # n_plots = length(heatmap_list)

    if (is.null(groups) || length(unique(groups)) == n_plots) {
        groups = 1:n_plots
    } else if (length(groups) == 1) {
        groups = rep(1, n_plots)
    } else if (length(groups) != n_plots) {
        stop("groups must have 1 value per plot")
    } else {
        groups = as.numeric(factor(groups, levels=unique(groups)))
    }

    n_groups = max(groups)
    # message(paste("groups:", paste(groups, collapse=', '), "\ntotal:", n_groups, "groups"))

    if (all(groups != 1:n_plots)) {
        for (i in 1:max(groups)) {
            group = which(groups == i)
            max_d = sum(vapply(heatmaps_list[group], function(x) max(x@value), integer(1)))
            for (index in group) {
                heatmaps_list[[index]]@max_value = max_d
            }
        }
    }

    # neater way of supplying group options
    dots = list(...)
    options[names(dots)] = dots
    opt_lengths = lengths(options)
    if (!all(opt_lengths %in% c(1, n_groups))) stop("supplied options must be length 1 or #groups")
    group_options = list()
    if (n_groups > 1) {
        for (i in 1:n_groups) {
            opt = options
            extra_opts = opt[lengths(opt) == n_groups]
            opt[names(extra_opts)] = lapply(extra_opts, function(x) unlist(x[i]))
            group_options[[i]] = opt
        }
    } else {
        group_options[[1]] = options
    }

    group_list = split(1:n_plots, groups)

    if (options$legend == TRUE) {
        widths = numeric(0)
        for (grp in group_list) {
            widths = c(widths, 0.2, rep(1, length(grp)))
        }
        mat = t(1:length(widths))
        widths = widths/sum(widths)
        layout(mat, widths)
    } else {
        layout(t(1:n_plots))
    }

    for (i in 1:n_groups) {
        go = group_options[[i]]
        grp = group_list[[i]]
        if(options$legend == TRUE) {
            message("plotting legend")
            par(mai=c(2, 2.8, 0.4, 0.1))
            plot_legend(heatmap_list[[grp[1]]]@max_value, go)

        }
        for (j in grp) {
            message("plotting heatmap")
            par(mai=c(2, 1, 0.4, 1))
            plotHeatmap(heatmap_list[[j]], go)
        }
    }
}


# additional functions

plot_legend <- function(max_value, options) {
        col_ramp <- colorRamp(.myColorPalette(options$color))
        ticks <- options$legend.ticks
        transf <- options$transform
        leg <- rep('', 256)
        leg[seq(1, 256, length.out=ticks)] <- formatC(seq(0, max_value, length.out=ticks), format='f', digits=2)
        plot(1, 1,
             type='n', bty='n',
             xaxt='n', yaxt='n',
             xlim=c(0,1), ylim=c(0,7),
             xaxs="i", yaxs="i",
             xlab='', ylab='')
        box(lwd = 6)
        color.legend(0, 0, 1, 7, legend=leg,
                     rect.col=rgb(col_ramp(transf(seq(0, 1, length.out=256)))/256),
                     align="lt", gradient='y', cex=options$cex.legend)
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
    return(colors.df[colorName,])
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

