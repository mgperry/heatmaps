#' Plot heatmaps for several patterns in DNA sequence
#'
#' @param seq DNAStringSet of equal width
#' @param patterns A vector or list of patterns
#' @param coords Heatmap coords
#' @param min.score Minimum score for PWM match
#' @param sigma Bandwith for smoothing kernel
#' @param output.ratio Output ratio of final image
#' @param output.size Output size of final image
#' @param options Heatmap plotting options
#' @param ... Additional Heatmap plotting options
#'
#' This function is a convenient wrapper for plotting many different patterns
#' for the same set of sequences. PatternHeatmap() is applied to the sequence
#' for each pattern in the list, they are passed to smooth() with the supplied
#' parameters and finally PlotHeatmapList().
#'
#' If fine-grained control is desired, or you want to mix other plot types, then
#' more information is available in the vignette.
#'
#' @seealso PatternHeatmap plotHeatmapList smooth
#' @export
#' @examples
#' data(HeatmapExamples)
#' plotPatternDensityMap(string_set, c("AT", "CG"), coords=c(-200, 200))
setGeneric(name="plotPatternDensityMap",
           def=function(seq, patterns, ...) standardGeneric("plotPatternDensityMap")
)

#' @describeIn plotPatternDensityMap Plot heatmaps for several patterns in DNA sequence
setMethod("plotPatternDensityMap",
signature(seq = "DNAStringSet"),
    function(seq, patterns, coords=NULL, min.score="80%", sigma=c(3, 3), output.ratio=c(1, 1), output.size=NULL, options=NULL, ...){

        if (is.null(options)) {
            options = heatmapOptions(...)
        }

        if (is.null(coords)) {
            coords = c(-width(seq[1])/2, width(seq[1])/2)
        }

        raw_hm = lapply(patterns, PatternHeatmap, seq=seq, min.score=min.score, coords=coords)
        if (!is.null(names(patterns))) {
            labels = ifelse(names(patterns) == "" | is.null(names(patterns)),
                            as.character(patterns),
                            names(patterns))
        } else {
            labels = as.character(patterns)
        }

        heatmaps <- list()
        for (i in 1:length(patterns)) {
            sm_hm = smooth(raw_hm[[i]], sigma=sigma, output.ratio=output.ratio, output.size=output.size, method="kernel")
            sm_hm@label = labels[i]
            heatmaps[[i]] = sm_hm
        }

        plotHeatmapList(heatmaps, groups=1, options)
    }
)

# setGeneric(
# name="plotPatternDensityMeta",
# def=function(seq, patterns, ...) {
#         standardGeneric("plotPatternDensityMeta")
#     }
# )

# setMethod("plotPatternDensityMeta",
# signature(seq = "DNAStringSet"),
# function(seq, patterns, coords=NULL, options=metaplotOptions(), ...){
#         if (is.null(coords)) {
#             coords = c(-width(seq[1])/2, width(seq[1])/2)
#         }

#         message("\nGetting oligonucleotide occurrence matrix...")
#         heatmaps = lapply(patterns, PatternHeatmap, seq=seq)
#         labels = sapply(patterns, function(x) if(class(x) == "PWM") x@name else x)

#         plotMetaRegion(heatmaps, groups=1, options)
#     }
# )

