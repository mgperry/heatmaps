setGeneric(
name="motifScanScores",
def=function(seq, PWM, ...){
    standardGeneric("motifScanScores")
    }
)

setMethod("motifScanScores",
signature(seq = "DNAStringSet", PWM = "matrix"),
function(seq, PWM, coords=NULL, label=NULL){

        if(!(length(unique(width(seq))) == 1)){
            stop("All sequences in the input DNAStringSet must have the same
            length!")
        }

        if (is.null(coords)) coords = c(0, length(seq[1]))

        if (is.null(label)) label=PWM@name

        scanning.score.list <- lapply(seq, function(x){
            PWMscoreStartingAt(pwm = PWM@matrix, subject = x,
                               starting.at = c(1:(length(seq[[1]]) - ncol(PWM@matrix) + 1)))
        })
        mat <- do.call(rbind, scanning.score.list)
        mat[mat < minScore(PWM)] <- minScore(PWM)

        # normalise to avoid < zero values
        max.score <- maxScore(PWM)
        min.score <- minScore(PWM)
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

