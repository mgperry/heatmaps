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
#' value, or more commonly as a percentage e.g. "80%". This is calculated based
#' on the difference between the minimum and maximum scores, not between zero and
#' the maximum score, since it is based on the log odds ratio. As of Bioconductor 3.5,
#' this behaviour is different to the `matchPWM` function in Biostrings.
#'
#' PatternHeatmaps often look much better after smoothing.
#'
#' @return A heatmap
#'
#' @seealso smoothHeatmap
#' @export
#' @examples
#' data(HeatmapExamples)
#' PatternHeatmap(string_set, "TA", coords=c(-100, 100), label="TA")
#' PatternHeatmap(string_set, tata_pwm, coords=c(-100, 100), min.score="80%", label="TATA PWM")
setGeneric(name="PatternHeatmap",
           def=function(seq, pattern, ...) standardGeneric("PatternHeatmap")
)

#' @describeIn PatternHeatmap Heatmap of sequence patterns from sequence and character
#' @importFrom Biostrings matchPattern startIndex DNAStringSet start maxScore minScore
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

        bp = width(seq[1])*length(seq)
        st = start(matchPattern(pattern, unlist(seq), fixed=FALSE))
        vec = numeric(bp)
        vec[st] = 1
        mat = matrix(vec, nrow=length(seq), byrow=TRUE)

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
        if (is.null(label)) label = "PWM"

        nc = nchar(min.score)
        if (substr(min.score, nc, nc) == "%") {
            percent = as.double(substr(min.score, 1L, nc-1L))
            min.score = pwm_score_percentile(pattern, percent/100)
        }

        bp = length(seq)*width(seq[1])
        st = start(matchPWM(subject=unlist(seq), pwm=pattern, min.score=min.score))
        vec = numeric(bp)
        vec[st] = 1
        mat = matrix(vec, nrow=length(seq), byrow=TRUE)

        hm = Heatmap(
            image=mat,
            scale=c(0,1),
            coords=as.integer(coords),
            nseq=length(seq),
            label=label)

        return(hm)
    }
)

pwm_score_percentile = function(pwm, percent) {
            min.value = minScore(pwm)
            max.value = maxScore(pwm)
            score.threshold = min.value + percent * (max.value-min.value)
}
