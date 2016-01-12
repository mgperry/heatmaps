setGeneric(
name="getPatternOccurrenceList",
def=function(regionsSeq, patterns, seqOrder = c(1:length(regionsSeq))){
        standardGeneric("getPatternOccurrenceList")
    }
)

# probably should work on only 1 pattern (collate later) and ignore seqOrder (handle at input)
setMethod("getPatternOccurrenceList",
signature(regionsSeq = "DNAStringSet"),
function(regionsSeq, patterns, seqOrder = c(1:length(regionsSeq))){
        
        if(!(length(unique(width(regionsSeq))) == 1)){
            stop("All sequences in the input DNAStringSet must have the same 
            length!")
        }
        if(!(length(seqOrder) == length(regionsSeq))){
            stop("The number of elements in 'seqOrder' must match the number 
            of input sequences in 'regionsSeq'!")
        }
        if(length(patterns) == 0){
            stop("At least one pattern needs to be specified!")
        }
        
        regionsSeq <- DNAStringSet(gsub("N", "+", regionsSeq))
        
        patterns.occurence.melted.list <- lapply(patterns, .get.pattern.occurence.melted.new, seq = regionsSeq, seqOrder = seqOrder)
        names(patterns.occurence.melted.list) <- patterns
        return(patterns.occurence.melted.list)
    }
)

.get.pattern.occurence.melted <- function(pattern, seq, seqOrder){
    pattern.matches <- vmatchPattern(pattern=pattern, subject=seq, fixed=FALSE)
    pattern.starts <- startIndex(pattern.matches)
    pattern.starts <- pattern.starts[seqOrder]
    non.zero.idx <- which(unlist(lapply(pattern.starts, length))>0)
    pattern.matrix.melt <- do.call(rbind, lapply(as.list(non.zero.idx),
        function(n) {
            data.frame(sequence = n, position = pattern.starts[[n]], value = 1)
        }))
    return(pattern.matrix.melt)
}

.get.pattern.occurence.melted.new <- function(pattern, seq, seqOrder){
    pattern.matches <- vmatchPattern(pattern=pattern, subject=seq, fixed=FALSE)
    pattern.starts <- startIndex(pattern.matches)
    pattern.starts <- pattern.starts[seqOrder]
    pattern.starts.ul <- unlist(pattern.starts)
    pattern.matrix.melt <- data.frame(
        sequence = rep(1:length(pattern.starts), lengths(pattern.starts)), 
        position = pattern.starts.ul, 
        value = 1)
    return(pattern.matrix.melt)
}





