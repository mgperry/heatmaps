#' Plot a list of heatmaps
#'
#' @param heatmap_list A list of Heatmaps
#' @param groups Optionally group heatmaps together
#' @param options Heatmap options
#' @param ... Additional options
#'
#' This function takes a list of one or more heatmaps and plots
#' them to a single image tiled horizontally.
#'
#' The "groups" argument specifies heatmaps to be grouped together and
#' plotted using the same display parameters and a unified scale. plotHeatmapList
#' will try to guess the best scale, either starting or finishing at zero, or symetrical
#' aronud zero - if this is not the desired behaviour, make sure the scales are
#' identical before the heatmaps are passed to the function.
#'
#' Options are specified as for plotHeatmap, but can be specified per group by passing
#' a list of options instead of a single vector. Note the difference between a length-2
#' character vector, c("Reds", "Blues"), and a list contatining two length-1 character
#' vectors: list("Reds", "Blues").
#'
#'
#' These are generally large, complex plots, so it can better to plot
#' straight to a file. PNG is preferred since pdf files generated can be
#' if the images are not downsized. The default settings are designed for plots
#' of about 10cm x 20cm per heatmap, but all of the relevant settigns can be tweaked
#' using the options. For display-quality images, it helps to increase the resolution
#' at to at least 150ppi, double the default of 72ppi on most systems.
#'
#' @return invisible(0)
#'
#' @seealso plotHeatmap heatmapOptions plot_legend
#' @importFrom graphics layout par
#' @export
#' @examples
#' data(HeatmapExamples)
#' plotHeatmapList(list(hm, hm2), groups=c(1,2), color=list("Reds", "Blues"))
plotHeatmapList = function(heatmap_list,
                           groups=1:length(heatmap_list),
                           options=heatmapOptions(),
                           ...) {

    # allow single heatmap argument
    if (class(heatmap_list) == "Heatmap") heatmap_list = list(heatmap_list)
    n_plots = length(heatmap_list)

    if (length(groups) == 1) {
        groups = rep(1, n_plots)
    } else {
        groups = as.numeric(factor(groups, levels=unique(groups)))
    }

    if (length(groups) != n_plots) stop("groups must have 1 value per plot")

    normalise_scales(heatmap_list, groups)

    dots = list(...)
    options[names(dots)] = dots

    options_list = expand_options(options, groups)

    params = get_device_params(options_list)
    layout(params$layout, params$width)
    par(cex=1) # can be changed by layout

    if (options$partition.legend) {
        par(mai=options_list[[1]]$legend.mai)
        plot_clusters(options)
    }

    for (i in 1:n_plots) {
        opts = options_list[[i]]
        hm = heatmap_list[[i]]
        message(paste("plotting heatmap", label(hm)))

        if(legend_left(opts)) {
            par(mai=opts$legend.mai)
            plot_legend(scale(hm), opts)
        }

        par(mai=opts$plot.mai)
        plotHeatmap(hm, opts)

        if(legend_right(opts)) {
            par(mai=opts$legend.mai[c(1,4,3,2)])
            plot_legend(scale(hm), opts)
        }
    }

    invisible(0)
}

legend_left = function(x) x$legend == TRUE && x$legend.pos == 'l'
legend_right = function(x) x$legend == TRUE && x$legend.pos == 'r'

normalise_scales = function(heatmaps, groups) {
    max = rep(-Inf, max(groups))
    min = rep(Inf, max(groups))
    for (i in seq_along(heatmaps)) {
        g = groups[i]
        s = scale(heatmaps[[i]])
        if (s[1] < min[g]) min[g] = s[1]
        if (s[2] > max[g]) max[g] = s[2]
    }
    for (i in seq_along(heatmaps)) {
        g = groups[i]
        scale(heatmaps[[i]]) = getScale(min[g], max[g])
    }
}

get_device_params = function(options_list) {
    if (options_list[[1]]$partition.legend) {
        width = options_list[[1]]$legend.width
    } else {
        width = numeric(0)
    }
    for (opts in options_list) {
        if (legend_left(opts)) width = c(width, opts$legend.width)
        width = c(width, 1)
        if (legend_right(opts)) width = c(width, opts$legend.width)
    }
    list(layout=t(1:length(width)), width=width)
}

expand_options = function(options, groups) {
    if (!all(lengths(Filter(is.list, options)) == max(groups))) {
        stop("Options supplied as lists must have length equal to the number of groups")
    }
    options_list = list()
    for (i in seq_along(groups)) {
        g = groups[i]
        opts = get_options(g, options)
        if (legend_left(opts)) {
            opts$legend = (min(which(groups == g)) == i) # plot only if first plot
        } else if (legend_right(opts)) {
            opts$legend = (max(which(groups == g)) == i) # plot  only if last plot
        }
        options_list[[i]] = opts
    }
    options_list
}

get_options = function(i, options) lapply(options, subset_if_list, i=i)

subset_if_list = function(i, x) { if (is.list(x)) x[[i]] else x }

#' Plot a color legend for a heatmap
#'
#' @param scale Numeric vector contain min and max for the scale
#' @param options heatmapOptions passed as a list
#'
#' This function plots a vertical color scale (or legend). With the default parameters,
#' it looks good at about 1/5 the width of a heatmap, about 1cm x 10cm. This
#' function only plots the legend, it does not set margin parameters.
#'
#' @return invisible(0)
#'
#' @seealso plotHeatmapList
#' @importFrom plotrix color.legend
#' @importFrom grDevices colorRamp rgb
#' @importFrom graphics box
#' @export
#' @examples
#' data(HeatmapExamples)
#' opts = heatmapOptions()
#' opts$color = "Rainbow"
#' par(mai=opts$legend.mai)
#' plot_legend(c(0,1), opts)
plot_legend <- function(scale, options) {
        if (length(options$color) == 1) {
            col_ramp = colorRamp(default_color(options$color)) # could check for function as input
        } else {
            col_ramp = colorRamp(options$color) # could check for function as input
        }
        ticks = options$legend.ticks
        values = seq(0, 1, length.out=256)
        if (!is.na(options$transform)) {
            values = options$transform(values)
        }
        color_palette = rgb(col_ramp(values)/256)
        align = ifelse(options$legend.pos=='l', 'lt', 'rb')
        leg = rep('', 256)
        leg[seq(1, 256, length.out=ticks)] = formatC(seq(scale[1], scale[2], length.out=ticks), format='g', digits=3)
        plot(1, 1,
             type='n', bty='n',
             xaxt='n', yaxt='n',
             xlim=c(0,1), ylim=c(0,7),
             xaxs="i", yaxs="i",
             xlab='', ylab='')
        box(lwd = options$box.width)
        color.legend(0, 0, 1, 7, legend=leg,
                     rect.col=color_palette,
                     align=align, gradient='y', cex=options$cex.legend)

        invisible(0)
}


#' Plot partition in a separate panel
#'
#' @param options heatmapOptions passed as a list
#'
#' Two heatmapOptions values are relevant:
#'
#' * partition Numeric vector containing relative sizes of the clusters
#' * colors Colors to use for clusters, additional colors are discarded
#'
#' This function plots a vertical color scale (or legend). With the default parameters,
#' it looks good at about 1/5 the width of a heatmap, about 1cm x 10cm. This
#' function only plots the legend, it does not set margin parameters.
#'
#' @return invisible(0)
#'
#' @seealso plotHeatmapList
#' @importFrom graphics box rect
#' @export
#' @examples
#' data(HeatmapExamples)
#' opts = heatmapOptions()
#' opts$partition = c(1,2,3,4)
#' par(mai=opts$legend.mai)
#' plot_clusters(opts)
plot_clusters <- function(options) {
    p = rev(options$partition)
    plot(1, 1,
         type='n', bty='n',
         xaxt='n', yaxt='n',
         xlim=c(0,1), ylim=c(0,sum(options$partition)),
         xaxs="i", yaxs="i",
         xlab='', ylab='')
    box(lwd = options$box.width)
    y_0 = 0
    for (i in seq_along(p)) {
        y_1 = y_0 + p[i]
        rect(0, y_0, 1, y_1, col=options$partition.col[i], lwd=options$box.width)
        y_0 = y_1
    }
    invisible(0)
}

