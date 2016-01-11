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



