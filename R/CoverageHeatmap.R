setGeneric(
name="CoverageHeatmap",
def=function(windows, track, ...){
        standardGeneric("CoverageHeatmap")
    }
)


setMethod("CoverageHeatmap",
signature(windows = "GenomicRanges", track="GenomicRanges"),
    function(windows, track, coords=NULL, weight=1, label=NULL, nbin=0) {
        # add bins, max.value ?

        if (is.null(label)) label=deparse(substitute(track))
        message(label)

        if (!(length(unique(width(windows))) == 1)) {
            stop("All sequences in the input DNAStringSet must have the same
            length!")
        }

        if (is.null(coords)) coords = c(0, width(windows[1]))

        cov = coverage(track, weight=weight)
        hm = CoverageHeatmap(windows, cov, coords, label, nbin)
        return(hm)
    }
)

setMethod("CoverageHeatmap",
signature(windows = "GenomicRanges", track="RleList"),
    function(windows, track, coords=NULL, label=NULL, nbin=0) {
        # add bins, max.value ?

        if (is.null(label)) label=deparse(substitute(track))

        if (!(length(unique(width(windows))) == 1)) {
            stop("All sequences in the input DNAStringSet must have the same
            length!")
        }

        if (is.null(coords)) coords = c(0, width(windows[1]))

        if (nbin==0) {
            list_of_rle = windowViews(windows, track)
            mat = do.call(rbind, lapply(list_of_rle, as.vector))
        } else if (nbin > 0) {
            sm = ScoreMatrixBin(track, windows, nbin)
            mat = as(sm, "matrix")
        }

        hm = new(
            "Heatmap",
            matrix=mat,
            scale=c(0,max(mat)),
            coords=as.integer(coords),
            nseq=length(windows),
            label=label)

        return(hm)
    }
)

windowViews = function(gr, cov) {
    ord = order(gr)
    rnk = rank(gr)
    gr = gr[ord]
    chrs = intersect(names(cov), as.character(seqlevels(gr)))
    myViews = Views(cov[chrs], as(gr, "RangesList")[chrs])
    scores = unlist(lapply(myViews, as.list), recursive=FALSE)
    scores = Map(function(x, s) if(s == "-") rev(x) else x, scores, as.vector(strand(gr)))
    return(scores[rnk])
}

