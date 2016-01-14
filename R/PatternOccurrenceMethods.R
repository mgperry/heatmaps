setGeneric(
name="getPatternOccurrence",
def=function(seq, pattern){
        standardGeneric("getPatternOccurrence")
    }
)

setMethod("getPatternOccurrence",
signature(seq = "DNAStringSet", pattern="character"),
    function(seq, pattern) {

        if (!(length(unique(width(seq))) == 1)) {
            stop("All sequences in the input DNAStringSet must have the same
            length!")
        }

        if (length(pattern) != 1) {
            stop("Only one pattern must to be specified!")
        }

        seq <- DNAStringSet(gsub("N", "+", seq)) # avoid matching Ns

        pattern.matches <- vmatchPattern(pattern=pattern, subject=seq, fixed=FALSE)
        pattern.starts <- startIndex(pattern.matches)
        pattern.starts.ul <- unlist(pattern.starts)
        sm = sparseMatrix(
            i = rep(1:length(pattern.starts), lengths(pattern.starts)),
            j = pattern.starts.ul,
            dims = c(length(seq), width(seq)[1]))
        return(sm)
    }
)

setMethod("getPatternOccurrence",
signature(seq = "DNAStringSet", pattern = "PWM"),
function(seq, pattern) {

        if (!(length(unique(width(seq))) == 1)) {
            stop("All sequences in the input DNAStringSet must have the same
            length!")
        }

        # Are Ns a problem?

        pattern.starts <- lapply(seq, function(x) {
            start(matchPWM(subject=x, pwm=pattern@matrix, min.score=pattern@min.score))
        })
        pattern.starts.ul <- unlist(pattern.starts)
        sm = sparseMatrix(
            i = rep(1:length(pattern.starts), lengths(pattern.starts)),
            j = pattern.starts.ul,
            dims = c(length(seq), width(seq)[1]))
        return(sm)
    }
)

