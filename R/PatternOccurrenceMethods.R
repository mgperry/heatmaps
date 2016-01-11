setGeneric(
name="getPatternOccurrenceList",
def=function(regionsSeq, patterns, seqOrder = c(1:length(regionsSeq)),
    useMulticore = FALSE, nrCores = NULL){
        standardGeneric("getPatternOccurrenceList")
    }
)

setMethod("getPatternOccurrenceList",
signature(regionsSeq = "DNAStringSet"),
function(regionsSeq, patterns, seqOrder = c(1:length(regionsSeq)),
    useMulticore = FALSE, nrCores = NULL){
        
        pt <- .Platform$OS.type
        if(useMulticore == TRUE){
            if(pt == "unix"){
                if("parallel" %in% rownames(installed.packages()) == FALSE){
                    stop("Cannot use multicore because package 'parallel' is 
                    not installed!")
                }else{
                    library(parallel)
                    if(is.null(nrCores)){
                        nrCores <- detectCores()
                    }
                }
            }else{
                useMulticore = FALSE
                warning("Multicore is not supported on non-Unix platforms! 
                Setting useMulticore=FALSE")
            }
        }
        
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
        
        if(useMulticore == TRUE){
            patterns.occurence.melted.list <- mclapply(as.list(patterns),
                function(x){
                    .get.pattern.occurence.melted(pattern = x, seq = regionsSeq,
                    seqOrder = seqOrder)
                }, mc.cores = nrCores)
        }else{
            patterns.occurence.melted.list <- lapply(as.list(patterns),
                function(x){
                    .get.pattern.occurence.melted(pattern = x, seq = regionsSeq,
                    seqOrder = seqOrder)
                })
        }
        
        names(patterns.occurence.melted.list) <- patterns
        return(patterns.occurence.melted.list)
    }
)
