setGeneric(
name="getPatternOccurrence",
def=function(seq, pattern, ...){
        standardGeneric("getPatternOccurrence")
    }
)

setMethod("getPatternOccurrence",
signature(seq = "DNAStringSet", pattern="character"),
    function(seq, pattern, coords=NULL) {

        if (!(length(unique(width(seq))) == 1)) {
            stop("All sequences in the input DNAStringSet must have the same
            length!")
        }

        if (length(pattern) != 1) {
            stop("Only one pattern must to be specified!")
        }

        if (is.null(coords)) coords = c(0, width(seq[1]))

        seq <- DNAStringSet(gsub("N", "+", seq)) # avoid matching Ns

        pattern.matches <- vmatchPattern(pattern=pattern, subject=seq, fixed=FALSE)
        pattern.starts <- startIndex(pattern.matches)
        pattern.starts.ul <- unlist(pattern.starts)

        sm = sparseMatrix(
            i = rep(1:length(pattern.starts), lengths(pattern.starts)),
            j = pattern.starts.ul,
            dims = c(length(seq), width(seq)[1]))

        mat = as.matrix(sm)*1

        hm = new(
            "Heatmap",
            matrix=mat,
            max_value=1,
            coords=as.integer(coords),
            nseq=length(seq),
            label=pattern)

        return(hm)
    }
)

setMethod("getPatternOccurrence",
signature(seq = "DNAStringSet", pattern = "PWM"),
function(seq, pattern, coords=NULL) {

        if (!(length(unique(width(seq))) == 1)) {
            stop("All sequences in the input DNAStringSet must have the same
            length!")
        }

        if (is.null(coords)) coords = c(0, width(seq[1]))

        pattern.starts <- lapply(seq, function(x) {
            start(matchPWM(subject=x, pwm=pattern@matrix, min.score=pattern@min.score))
        })
        pattern.starts.ul <- unlist(pattern.starts)

        sm = sparseMatrix(
            i = rep(1:length(pattern.starts), lengths(pattern.starts)),
            j = pattern.starts.ul,
            dims = c(length(seq), width(seq)[1]))

        mat = as.matrix(sm)*1

        hm = new(
            "Heatmap",
            matrix=sm,
            max_value=mat,
            coords=as.integer(coords),
            nseq=length(seq),
            label=pattern)
        return(hm)
    }
)

setGeneric(
name="CoverageHeatmap",
def=function(windows, track, ...){
        standardGeneric("CoverageHeatmap")
    }
)


setMethod("CoverageHeatmap",
signature(windows = "GenomicRanges", track="GenomicRanges"),
    function(windows, track, coords=NULL, weight=1, label=NULL) {
        # add bins, max.value ?

        if (is.null(label)) label=deparse(substitute(track))
        message(label)

        if (!(length(unique(width(windows))) == 1)) {
            stop("All sequences in the input DNAStringSet must have the same
            length!")
        }

        if (is.null(coords)) coords = c(0, width(windows[1]))

        cov = coverage(track, weight=weight)
        list_of_rle = windowViews(windows, cov)

        mat = do.call(rbind, lapply(list_of_rle, as.vector))

        hm = new(
            "Heatmap",
            matrix=mat,
            max_value=max(mat),
            coords=as.integer(coords),
            nseq=length(windows),
            label=label)
        return(hm)
    }
)

setMethod("CoverageHeatmap",
signature(windows = "GenomicRanges", track="RleList"),
    function(windows, track, coords=NULL, label=NULL) {
        # add bins, max.value ?

        if (is.null(label)) label=deparse(substitute(track))

        if (!(length(unique(width(windows))) == 1)) {
            stop("All sequences in the input DNAStringSet must have the same
            length!")
        }

        if (is.null(coords)) coords = c(0, width(windows[1]))

        list_of_rle = windowViews(windows, track)

        mat = do.call(rbind, lapply(as.vector(list_of_rle)))

        hm = new(
            "Heatmap",
            matrix=mat,
            max_value=max(mat),
            coords=as.integer(coords),
            nseq=length(windows),
            label=label)
        return(hm)
    }
)

windowViews = function(gr, coverage) {
    ord = order(gr)
    rev_ord = seq(1, length(gr))[ord]
    gr = gr[ord]
    chrs = intersect(names(coverage), as.character(seqlevels(gr)))
    myViews = Views(coverage[chrs], as(gr, "RangesList")[chrs])
    scores = unlist(lapply(myViews, as.list), recursive=FALSE)
    return(scores[rev_ord])
}

