setGeneric(
name="PWMScanHeatmap",
def=function(seq, pwm, ...){
    standardGeneric("motifScanScores")
    }
)

setMethod("PWMScanHeatmap",
signature(seq = "DNAStringSet", pwm = "matrix"),
function(seq, pwm, coords=NULL, label=NULL){

        if(!(length(unique(width(seq))) == 1)){
            stop("All sequences in the input DNAStringSet must have the same
            length!")
        }

        if (is.null(coords)) coords = c(0, length(seq[1]))

        scanning.score.list <- lapply(seq, function(x){
           pwmscoreStartingAt(pwm = pwm, subject = x,
                               starting.at = c(1:(length(seq[[1]]) - ncol(pwm + 1))))
        })
        mat <- do.call(rbind, scanning.score.list)
        mat[mat < minScore(pwm)] <- minScore(pwm)

        # normalise to avoid < zero values
        max.score <- maxScore(pwm)
        min.score <- minScore(pwm)
        mat = (mat-min.score)/(max.score-min.score)*100

        hm = new(
            "Heatmap",
            matrix=mat,
            max_value=max(mat),
            coords=as.integer(coords),
            nseq=length(seq),
            label=label)

        return(hm)
    }
)

