setGeneric(
name="PatternHeatmap",
def=function(seq, pattern, ...){
        standardGeneric("PatternHeatmap")
    }
)

setMethod("PatternHeatmap",
signature(seq = "DNAStringSet", pattern="character"),
    function(seq, pattern, coords=NULL, min.score=NULL, label=NULL) {

        if (!(length(unique(width(seq))) == 1)) {
            stop("All sequences in the input DNAStringSet must have the same
            length!")
        }

        if (length(pattern) != 1) {
            stop("Only one pattern must to be specified!")
        }

        if (is.null(coords)) coords = c(0, width(seq[1]))
        if (is.null(label)) label = pattern

        seq <- DNAStringSet(gsub("N", "+", seq)) # avoid matching Ns

        pattern.matches <- vmatchPattern(pattern=pattern, subject=seq, fixed=FALSE)
        pattern.starts <- startIndex(pattern.matches)
        pattern.starts.ul <- unlist(pattern.starts) + floor(length(pattern)/2)

        sm = sparseMatrix(
            i = rep(1:length(pattern.starts), lengths(pattern.starts)),
            j = pattern.starts.ul,
            dims = c(length(seq), width(seq)[1]))

        mat = as.matrix(sm)*1

        hm = new(
            "Heatmap",
            matrix=mat,
            scale=c(0,1),
            coords=as.integer(coords),
            nseq=length(seq),
            label=label)

        return(hm)
    }
)

setMethod("PatternHeatmap",
signature(seq = "DNAStringSet", pattern = "matrix"),
function(seq, pattern, coords=NULL, min.score="80%", label=NULL) {

        if (!(length(unique(width(seq))) == 1)) {
            stop("All sequences in the input DNAStringSet must have the same
            length!")
        }

        if (is.null(coords)) coords = c(0, width(seq[1]))
        if (is.null(label)) label = "pwm"

        pattern.starts <- lapply(seq, function(x) {
            start(matchPWM(subject=x, pwm=pattern, min.score=min.score))
        })
        pattern.starts.ul <- unlist(pattern.starts) + floor(ncol(pattern)/2)

        sm = sparseMatrix(
            i = rep(1:length(pattern.starts), lengths(pattern.starts)),
            j = pattern.starts.ul,
            dims = c(length(seq), width(seq)[1]))

        mat = as.matrix(sm)*1

        hm = new(
            "Heatmap",
            matrix=mat,
            scale=c(0,max(mat)),
            coords=as.integer(coords),
            nseq=length(seq),
            label=label)
        return(hm)
    }
)

