#' Plot a Heatmap object to the device
#'
#' @param heatmap A heatmap object
#' @param options A list containing plotting options
#' @param ... Used for passing individual options
#'
#' This function will take a heatmap and plot it to the device with
#' the specified options. Options can be passed together in a list
#' or individually as additional arguments. If passing options as a list,
#' it's best to first create a list containing the default settings using
#' heatmapOptions() andmethod  then setting options individually.
#'
#' plotHeatmap() does not control device settings at all, these can be
#' set using plotHeatmapList() and the relevant options in heatmapOptions()
#'
#' See ?heatmapOptions for a full list of options.
#'
#' @return invisible(0)
#'
#' @export
#' @seealso heatmapOptions plotHeatmapList
#' @importFrom grDevices colorRamp colorRampPalette
#' @importFrom graphics abline axis box grconvertX grconvertY lines text
#' @examples
#'
#' data(HeatmapExamples)
#' plotHeatmap(hm, color="Blues")
setGeneric(name="plotHeatmap",
           def=function(heatmap, options=NULL, ...) standardGeneric("plotHeatmap")
)

#' @describeIn plotHeatmap Plot a Heatmap object to the device
#' @export
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
        val = options$transform(image(heatmap))
        scale_t = options$transform(scale(heatmap))
        breaks = seq(scale_t[1], scale_t[2], length.out=257)
    } else {
        val = image(heatmap)
        breaks = seq(scale(heatmap)[1], scale(heatmap)[2], length.out=257)
    }

    val[val < scale(heatmap)[1]] = scale(heatmap)[1]
    val[val > scale(heatmap)[2]] = scale(heatmap)[2]

    image(xm(heatmap), ym(heatmap), z=t(val)[,nrow(val):1],
          col=col_palette(256), breaks=breaks,
          xlim=c(0.5,width(heatmap)-0.5), # remove fluffy edges
          xlab="", ylab="",
          axes=FALSE,
          useRaster=TRUE)

    if (options$box.width > 0) box(lwd = options$box.width)

    if(options$x.ticks == TRUE) {
        x.ticks = make_x_ticks(heatmap@coords, options$x.tick.labels)
        axis(1, at=x.ticks, labels=names(x.ticks), cex.axis=options$cex.axis, lwd=options$box.width, tcl=options$tcl, padj=options$padj)
    }

    if(length(options$partition) > 1 && options$partition.lines == TRUE) {
        line_coords = heatmap@nseq * (1 - cumsum(options$partition/sum(options$partition)))
        abline(h=line_coords[1:(length(line_coords) - 1)], lwd=options$box.width)
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
             labels=paste(round(scale.length), options$scale.label, sep=''),
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

    if(!is.na(options$hook)){
        options$hook()
    }

    invisible(0)
})

#' Generate default options for a Heatmap
#'
#' @param ... options to set manually
#'
#' Guide to Heatmap options
#'
#' This is an reference to all the possible options for plotting heatmaps.
#' Some options are handled by heatmaps functions (either plotHeatmap or plotHeatmapList),
#' others are passed directly to plotting functions. Further explanation is available in the
#' vignette. Arguments are numeric if not otherwise stated.
#'
#' color: A vector of colors or a default color, see ?default_color. plotHeatmap will
#'        interpolate between these colors to form a scale.
#'
#' box.width: width of box around the heatmap, passed to box()
#'
#' x.ticks: Logical, plot x axis ticks
#'
#' x.tick.labels: Character, labels to use for x ticks, (default blank)
#'
#' tcl: Length of x axis ticks
#'
#' padj: Vertical adjustment of x axis labels
#'
#' cex.axis: cex for axis labels
#'
#' scale: Logical, Plot scale or not
#'
#' scale.label: Character, label for scale
#'
#' scale.lwd: Width for line around scale
#'
#' cex.scale: Cex for Scale
#'
#' label: Logical, plot label or not
#'
#' label.xpos: x position for label, from left
#'
#' label.ypos: y position for label, from top
#'
#' cex.label: cex for axis labels
#'
#' label.col: Color for label, white is often useful for dark plots
#'
#' legend: Logical, plot legend (scale indicating values for colors)
#'
#' legend: Color for label, white is often useful for dark plots
#'
#' legend.pos: Character, position of legend relative to heatmap: 'l' for left, 'r' for right
#'
#' legend.ticks: Number of ticks to use on legend.
#'
#' cex.legend: cex to use for legend marks
#'
#' refline: Logical, Draw dashed line at coords = 0
#' label: Logical, plot label or not
#'
#' label.xpos: x position for label, from left
#'
#' label.ypos: y position for label, from top
#'
#' cex.label: cex for axis labels
#'
#' label.col: Color for label, white is often useful for dark plots
#'
#' legend: Logical, plot legend (scale indicating values for colors)
#'
#' legend: Color for label, white is often useful for dark plots
#'
#' legend.pos: Character, position of legend relative to heatmap: 'l' for left, 'r' for right
#'
#' legend.ticks: Number of ticks to use on legend.
#'
#' cex.legend: cex to use for legend marks
#'
#' refline: Logical, Draw dashed line at coords = 0
#'
#' refline.width: Width of reference line
#'
#' transform: Function to transform values before plotting
#'
#' plot.mai: Length-4 numeric, margins around plot
#'
#' legend.mai: Length-4 numeric, margins around legend
#'
#' partition: Numeric, relative sizes of clusters
#'
#' partition.lines: Logical, plot lines delineating clusters
#'
#' partition.legend: Logical, plot cluster legend in HeatmapList
#'
#' partition.col: Character, colours to use for plotting clusters. Defaults to RColorBrewer's Set1
#'
#' hook: Function called after plotting is complete.
#'
#' @return a list containing the specified options
#' @importFrom RColorBrewer brewer.pal
#' @seealso plotHeatmap plotHeatmapList
#' @export
#' @examples
#'
#' myOptions = heatmapOptions()
#' myOptions$color = "Reds"
#' # plotHeatmap(hm, options=myOptions)
heatmapOptions = function(...) {
    opts = list(...)
    default = list(
        color='Blues',
        box.width=2.5,
        x.ticks=TRUE,
        x.tick.labels='',
        tcl=-0.5,
        padj=0,
        cex.axis=1.5,
        scale=FALSE,
        scale.label='bp',
        scale.lwd=6,
        cex.scale=1.5,
        label=TRUE,
        label.xpos=0.1,
        label.ypos=0.1,
        cex.label=2.5,
        label.col='black',
        legend=FALSE,
        legend.width=0.2,
        legend.pos='l',
        legend.ticks=5,
        cex.legend=1.2,
        refline=TRUE,
        refline.width=1,
        transform=NA,
        plot.mai=c(0.6, 0.3, 0.1, 0.3),
        legend.mai=c(0.6, 0.6, 0.1, 0.05),
        partition=1,
        partition.lines=FALSE,
        partition.legend=FALSE,
        partition.col=brewer.pal(9, "Set1"),
        hook=NA
    )
    bad_names = names(opts)[!names(opts) %in% names(default)]
    if (length(bad_names) > 0) warning(paste("Arguments", bad_names, "are not heatmap options"))
    default[names(opts)] = opts
    return(default)
}

make_x_ticks = function(coord, label) {
    if (coord[1] < 0) {
        labs = c(1*coord[1], "", 0, "", coord[2])
        labs[c(1,3,5)] = paste0(labs[c(1,3,5)], label)
        x_ticks= cumsum(c(0.5, -coord[1]/2, -coord[1]/2, coord[2]/2-1, coord[2]/2))
    } else {
        labs = c(coord[1], "", coord[2])
        labs[c(1,3)] = paste0(labs[c(1,3)], label)
        x_ticks= cumsum(c(0.5, ((coord[2]-coord[1])/2)-0.5, ((coord[2]-coord[1])/2)-0.5))
    }
    names(x_ticks) = labs
    x_ticks
}

color_palettes = list(
  Rainbow=rev(rainbow(9, start=0, end=4/6)),
  Rasta=c("#E20D0D", "#F0D817", "#3EB308")
)

#' Predifined color palettes from RColorBrewer + Rainbow
#'
#' @param col Character, RColorBrewer colorscheme or "Rainbow"
#'
#' This function provides a convenient function to all color palettes from
#' RColorBrewer, and a better version of R's rainbow function (specifically
#' rev(rainbow(9, start=0, end=4/6)), so it starts blue and ends with red).
#'
#' @return character, a length-9 color palette
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices rainbow
#' @export
#' @examples
#' default_color("Blues")
#' default_color("Rainbow")
default_color = function(col) {
    if (col %in% names(color_palettes)) {
        palette = color_palettes[[col]]
    } else {
       palette = tryCatch(brewer.pal(9, col), error=color_error)
    }
    palette
}

color_error = function(e) "Incorrect color specification: refer to ?default_color"


