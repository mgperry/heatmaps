#' Generate a Heatmap of coverage
#'
#' @param windows A set of GRanges of equal length
#' @param track A GRanges or RleList object specifying coverage
#' @param coords Co-ordinates for the heatmap, defaults to c(0, width(windows))
#' @param weight Passed to coverage(track) constructor if class(track) == "GRanges"
#' @param label Label for the heatmap
#' @param nbin If set, number of bins to use across each window
#' @param ... additional arguments used by methods
#'
#' This function generates a Heatmap object from a set of windows and an
#' object containing genome-wide information about coverage. Either a GRanges
#' or an RleList can be used. In the former case, the "weight" paramter is
#' passed directly to the `coverage` function. If nbin is set, binned coverage
#' is calculated which will save memory and time when plotting and average out
#' varible data.
#'
#' If the coverage track contains negative values, then the scale will be
#' centered on zero, ie. c(-max(abs(image(hm))), max(abs(image(hm)))). This
#' makes more sense for most color schemes which are centered on zero, and
#' avoids misleading plots where either positive or negative values are
#' over-emphasised. See ?getScale for details. The scale can be manually reset
#' if desired using the "scale" method.
#'
#' @return A Heatmap object
#'
#' @importFrom GenomicRanges coverage strand
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom IRanges Views revElements tile mean
#' @export
#' @examples
#' data(HeatmapExamples)
#' CoverageHeatmap(windows, rle_list, coords=c(-100, 100), label="Example")
setGeneric(name="CoverageHeatmap",
           def=function(windows, track, ...) standardGeneric("CoverageHeatmap")
)


#' @describeIn CoverageHeatmap Heatmap of Coverage from 2 GRanges
#' @export
setMethod("CoverageHeatmap",
signature(windows = "GenomicRanges", track="GenomicRanges"),
    function(windows, track, coords=NULL, weight=1, label=NULL, nbin=0) {

        w = unique(width(windows))
        if (length(w) != 1) stop("All ranges in windows must have the same length!")

        if (is.null(coords)) coords = c(0, width(windows[1]))
        if (is.null(label)) label=deparse(substitute(track))

        cov = coverage(track, weight=weight)
        hm = CoverageHeatmap(windows, cov, coords, label, nbin)
        return(hm)
    }
)

#' @describeIn CoverageHeatmap Heatmap of Coverage from GRanges + RleList
#' @export
setMethod("CoverageHeatmap",
signature(windows = "GenomicRanges", track="RleList"),
    function(windows, track, coords=NULL, label=NULL, nbin=0) {

        w = unique(width(windows))
        if (length(w) != 1) stop("All ranges in windows must have the same length!")

        if (is.null(label)) label=deparse(substitute(track))
        if (is.null(coords)) coords = c(0, width(windows[1]))

        neg = strand(windows) == "-"

        if (nbin==0) {
            rle_list = track[windows]
            rle_list[neg] = revElements(rle_list[neg])
            mat = as.matrix(rle_list)
        } else if (nbin > 0) {
            tiles = tile(windows, nbin)
            tiles[neg] = revElements(tiles[neg]) # window direction does not matter here
            mat = matrix(mean(track[unlist(tiles)]), ncol=nbin, byrow=TRUE)
        }

        hm = Heatmap(
            image=mat,
            scale=getScale(min(mat), max(mat)),
            coords=as.integer(coords),
            nseq=length(windows),
            label=label)

        return(hm)
    }
)

