#' Generate a Heatmap of PWM Scores in DNA sequnce
#'
#' @param seq A DNAString of equal length
#' @param pwm A PWM
#' @param coords Co-ordinates for the heatmap, defaults to c(0, width(windows))
#' @param label Label for the heatmap
#' @param ... additional arguments used by methods
#'
#' This function creates a heatmap where each point is the score of a PWM match
#' starting from that position, which can visualise regions of enrichment or exclusion
#' of certain motifs
#'
#' @return A heatmap
#'
#' @seealso PatternHeatmap
#' @export
#' @examples
#' data(HeatmapExamples)
#' PatternHeatmap(string_set, tata_pwm, coords=c(-100, 100), label="TATA Scan")
setGeneric(name="PWMScanHeatmap",
           def=function(seq, pwm, ...) standardGeneric("PWMScanHeatmap")
)

#' @describeIn PWMScanHeatmap Heatmap of PWM Scores
#' @importFrom Biostrings PWMscoreStartingAt minScore maxScore
#' @export
setMethod("PWMScanHeatmap",
signature(seq = "DNAStringSet", pwm = "matrix"),
function(seq, pwm, coords=NULL, label=""){

        if(!(length(unique(width(seq))) == 1)){
            stop("All sequences in the input DNAStringSet must have the same
            length!")
        }

        if (is.null(coords)) coords = c(0, length(seq[[1]]))

        # new fast way
        w = width(seq[1])
        columns = w - ncol(pwm)
        bp = length(seq)*width(seq)[[1]]
        scores = PWMscoreStartingAt(pwm, unlist(seq), starting.at = c(1:(bp-ncol(pwm))))
        mat = matrix(c(scores, rep(0, ncol(pwm))), ncol=width(seq), byrow=TRUE)

        # normalise to avoid < zero values
        max.score <- maxScore(pwm)
        min.score <- minScore(pwm)
        mat = (mat-min.score)/(max.score-min.score)*100
        mat[,columns] = 50 # set irrelevant columns to average

        hm = Heatmap(
            image=mat,
            scale=c(0, 100),
            coords=as.integer(coords),
            nseq=length(seq),
            label=label)

        return(hm)
    }
)

