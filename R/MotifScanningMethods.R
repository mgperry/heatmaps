setGeneric(
name="motifScanHits",
def=function(regionsSeq, motifPWM, minScore = "80%",
    seqOrder = c(1:length(regionsSeq))){
        standardGeneric("motifScanHits")
    }
)

setMethod("motifScanHits",
signature(regionsSeq = "DNAStringSet", motifPWM = "matrix"),
function(regionsSeq, motifPWM, minScore = "80%",
    seqOrder = c(1:length(regionsSeq))){
    
        if(!(length(unique(width(regionsSeq))) == 1)){
            stop("All sequences in the input DNAStringSet must have the same 
            length!")
        }
        if(!(length(seqOrder) == length(regionsSeq))){
            stop("The number of elements in 'seqOrder' must match the number 
            of input sequences in 'regionsSeq'!")
        }

        pwm.match <- lapply(regionsSeq[seqOrder], function(x){
            .scan.sequence.with.pwm(pwm = motifPWM, seq = x, minScore=minScore)
        })
        
        data.frame(sequence = rep(c(1:length(pwm.match)),
        times = unlist(lapply(pwm.match, length))),
        position = unlist(pwm.match), value = rep(1, times = length(unlist(pwm.match))))
    
    }
)


setGeneric(
name="motifScanScores",
def=function(regionsSeq, motifPWM, seqOrder = c(1:length(regionsSeq)),
    asPercentage = TRUE){
    standardGeneric("motifScanScores")
    }
)

setMethod("motifScanScores",
signature(regionsSeq = "DNAStringSet", motifPWM = "matrix"),
function(regionsSeq, motifPWM, seqOrder = c(1:length(regionsSeq)),
    asPercentage = TRUE){
        
        if(!(length(unique(width(regionsSeq))) == 1)){
            stop("All sequences in the input DNAStringSet must have the same 
            length!")
        }
        if(!(length(seqOrder) == length(regionsSeq))){
            stop("The number of elements in 'seqOrder' must match the number 
            of input sequences in 'regionsSeq'!")
        }
        
        regionsSeq <- regionsSeq[seqOrder]
        scanning.score.list <- lapply(regionsSeq, function(x){
            PWMscoreStartingAt(pwm = motifPWM, subject = x, starting.at =
            c(1:(length(regionsSeq[[1]]) - ncol(motifPWM) + 1)))
        })
        scanning.score.matrix <- do.call(rbind, scanning.score.list)
        scanning.score.matrix[scanning.score.matrix < minScore(motifPWM)] <- minScore(motifPWM)

        if(asPercentage){
            max.score <- maxScore(motifPWM)
            min.score <- minScore(motifPWM)
            return((scanning.score.matrix-min.score)/(max.score-min.score)*100)
        }else{
            return(scanning.score.matrix)
        }
    
    }
)

