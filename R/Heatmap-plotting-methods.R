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

    xm <- 0:ncol(heatmap@matrix)
    ym <- 0:nrow(heatmap@matrix)
    val <- options$transform(heatmap@matrix)
    val[val < 0] = 0
    val[val > heatmap@max_value] = heatmap@max_value

    image(xm, ym, z=t(val),
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
    usr = list(...)
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
        transform=function(x) x)
    def[names(usr)] = usr
    return(def)
}

plotHeatmapList = function(heatmap_list, groups=NULL, options=heatmapOptions(), ...) {
    if (class(heatmap_list) == "Heatmap") heatmap_list = list(heatmap_list) # allow single heatmap argument
    n_plots = length(heatmap_list)

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
        for (i in 1:length(group_list)) {
            if (group_options[[i]]$legend.pos == 'l') {
                widths = c(widths, 0.2, rep(1, length(group_list[[i]])))
            } else if (group_options[[i]]$legend.pos == 'r') {
                widths = c(widths, rep(1, length(group_list[[i]])), 0.2)
            }
        }
        mat = t(1:length(widths))
        widths = widths/sum(widths)
        layout(mat, widths)
    } else {
        layout(t(1:n_plots))
    }

    par(cex=1) # can be changed by layout

    for (i in 1:n_groups) {
        go = group_options[[i]]
        grp = group_list[[i]]
        if(go$legend == TRUE && go$legend.pos == 'l') {
            message("plotting legend")
            par(mai=c(2, 2.8, 0.4, 0.1))
            plot_legend(heatmap_list[[grp[1]]]@max_value, go)
        }
        for (j in grp) {
            message("plotting heatmap")
            par(mai=c(2, 1.2, 0.4, 1.2))
            plotHeatmap(heatmap_list[[j]], go)
        }
        if(go$legend == TRUE && go$legend.pos == 'r') {
            message("plotting legend")
            par(mai=c(2, 2.8, 0.4, 0.1))
            plot_legend(heatmap_list[[grp[1]]]@max_value, go)
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

