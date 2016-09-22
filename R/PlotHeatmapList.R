#' PLot a list of Heatmaps
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
#' Options are specified as for plotHeatmap, but can be specified per group.
#' Some parameters can be passed as vectors, but care must be taken for parameters
#' which are themselves vectors of length > 1, for example color palettes, which must
#' be specified as lists. Examples of this are documented in the vignette.
#'
#' These are generally large, complex plots, so it can better to plot
#' straight to a file. PNG is preferred since pdf files generated can be
#' if the images are not downsized. The default settings are designed for plots
#' of about 5cm x 10cm per heatmap, but all of the relevant settigns can be tweaked
#' using the options. For display-quality images, it helps to increase the resolution
#' at to at least 150ppi, double the default of 72ppi on most systems.
#'
#' @seealso plotHeatmap heatmapOptions plot_legend
#' @export
#' @examples
#' data(HeatmapExamples)
#' plotHeatmapList(list(hm, hm2), groups=c(1,2), color=list("Reds", "Blues"))
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

    for (i in 1:max(groups)) {
        group = which(groups == i)
        # test if scales are equal
        scales = vapply(heatmap_list[group], function(x) scale(x), numeric(2))
        if (!(length(unique(scales[1,])) == 1) && (length(unique(scales[2,])) == 1)) {
            for (index in group) {
                scale(heatmap_list[[index]]) = get_scale(min(scales), max(scales))
            }
        }
    }

    # neater way of supplying group options
    dots = list(...)
    options[names(dots)] = dots
    opt_lengths = lengths(options)
    if (!all(opt_lengths %in% c(1, n_groups))) stop("supplied options must be length 1 or #groups")
    group_options = list()
    for (i in 1:n_groups) {
        opt = options
        extra_opts = opt[lengths(opt) == n_groups]
        for (n in names(extra_opts)) {
            opt[[n]] = extra_opts[[n]][i]
        }
        group_options[[i]] = lapply(opt, safe_unlist)
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
            par(mai=go$legend.mai[[1]])
            plot_legend(scale(heatmap_list[[grp[1]]]), go)
        }
        for (j in grp) {
            message(paste("plotting heatmap", go$label))
            par(mai=go$plot.mai[[1]])
            plotHeatmap(heatmap_list[[j]], go)
        }
        if(go$legend == TRUE && go$legend.pos == 'r') {
            par(mai=go$legend.mai[[1]])
            plot_legend(scale(heatmap_list[[grp[1]]]), go)
        }
    }
}

safe_unlist = function(x) {
    for (i in seq_along(x)) {
        if (is.list(x[[i]])) {
            x[[i]] = x[[i]][[1]]
        }
    }
    x
}

#' Plot a color legend for a heatmap
#'
#' @param scale Numeric vector contain min and max for the scale
#' @param options heatmapOptions passed as a list
#'
#' This function plots a vertical color scale (or legend). With the default parameters,
#' it looks good at about 1/5 the width of a heatmap, about 1cm x 10cm. This
#' function only plots the legend, it does not set margin parameters.
#'
#' @seealso plotHeatmapList
#' @importFrom plotrix color.legend
#' @examples
#' data(HeatmapExamples)
#' opts = heatmapOptions
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
        leg[seq(1, 256, length.out=ticks)] = formatC(seq(scale[1], scale[2], length.out=ticks), format='f', digits=2)
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
}

