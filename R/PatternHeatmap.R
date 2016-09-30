#' Generate a Heatmap of patterns in DNA sequence
#'
#' @param seq A DNAString of equal length
#' @param pattern A nucleotide pattern or PWM
#' @param coords Co-ordinates for the heatmap, defaults to c(0, width(windows))
#' @param min.score Minimum score for PWM match
#' @param label Label for the heatmap
#' @param ... additional arguments used by methods
#'
#' This function creates a Heatmap from a set of DNA sequences. The resulting
#' heatmap will be binary, with 1 representing a match and 0 otherwise. Patterns
#' can be specified as a character vectore, eg. "CTCCC", or as a PWM. These
#' arguments are passed to Biostrings functions, `vmatchPattern` and `matchPWM`.
#' Character arguments can contain standard ambiguity codes. PWMs must be 4 by n
#' matricies with columns names ACGT. "min.score" is specified either as an absolute
#' value, or more commonly as a percentage e.g. "80%". Refer to Biostrings documentation
#' for details.
#'
#' PatternHeatmaps often look much better after smoothing.
#'
#' @seealso smooth
#' @export
#' @examples
#' data(HeatmapExamples)
#' PatternHeatmap(string_set, "TA", coords=c(-100, 100), label="TA")
#' PatternHeatmap(string_set, tata_pwm, coords=c(-100, 100), min.score="80%", label="TATA PWM")
setGeneric(name="PatternHeatmap",
           def=function(seq, pattern, ...) standardGeneric("PatternHeatmap")
)

#' @describeIn PatternHeatmap Heatmap of sequence patterns from sequence and character
#' @importFrom Biostrings vmatchPattern startIndex DNAStringSet start
#' @importFrom Matrix sparseMatrix
#' @export
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

        hm = Heatmap(
            image=mat,
            scale=c(0,1),
            coords=as.integer(coords),
            nseq=length(seq),
            label=label)

        return(hm)
    }
)

#' @describeIn PatternHeatmap Heatmap of sequence patterns from sequence and matrix
#' @importFrom Biostrings matchPWM
#' @export
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

        hm = Heatmap(
            image=mat,
            scale=c(0,1),
            coords=as.integer(coords),
            nseq=length(seq),
            label=label)
        return(hm)
    }
)

