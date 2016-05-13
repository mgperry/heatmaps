setGeneric(
    name="plotPatternDensityMap",
    def=function(seq, patterns, ...) {
            standardGeneric("plotPatternDensityMap")
    }
)

setMethod("plotPatternDensityMap",
signature(seq = "DNAStringSet"),
    function(seq, patterns, coords=NULL, min.score="80%", nbin=NULL, vsmooth=NULL, hsmooth=NULL, bw=NULL, options=NULL, ...){

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
            sm_hm = smooth(raw_hm[[i]], nbin, vsmooth, hsmooth, bw)
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

