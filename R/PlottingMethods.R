setGeneric(
name="plotPatternDensityMap",
def=function(regionsSeq, patterns, seqOrder = c(1:length(regionsSeq)),
    flankUp = NULL, flankDown = NULL, nBin = NULL, bandWidth = NULL,
    color = "blue", transf = NULL, xTicks = NULL, xTicksAt = NULL, xLabel = "",
    yTicks = NULL, yTicksAt = NULL, yLabel = "", cexAxis = 8, plotScale = TRUE,
    scaleLength = NULL, scaleWidth = 15, addPatternLabel = TRUE, cexLabel = 8,
    labelCol = "black", addReferenceLine = TRUE, plotColorLegend = TRUE,
    outFile = "PatternDensityMap", plotWidth = 2000, plotHeight = 2000,
    useMulticore = FALSE, nrCores = NULL){
        standardGeneric("plotPatternDensityMap")
    }
)

setMethod("plotPatternDensityMap",
signature(regionsSeq = "DNAStringSet"),
function(regionsSeq, patterns, seqOrder = c(1:length(regionsSeq)),
    flankUp = NULL, flankDown = NULL, nBin = NULL, bandWidth = NULL,
    color = "blue", transf = NULL, xTicks = NULL, xTicksAt = NULL, xLabel = "",
    yTicks = NULL, yTicksAt = NULL, yLabel = "", cexAxis = 8, plotScale = TRUE,
    scaleLength = NULL, scaleWidth = 15, addPatternLabel = TRUE, cexLabel = 8,
    labelCol = "black", addReferenceLine = TRUE, plotColorLegend = TRUE,
    outFile = "PatternDensityMap", plotWidth = 2000, plotHeight = 2000,
    useMulticore = FALSE, nrCores = NULL){
        
        pt <- .Platform$OS.type
        if(useMulticore == TRUE){
            if(pt == "unix"){
                if("parallel" %in% rownames(installed.packages()) == FALSE){
                    stop("Cannot use multicore because package 'parallel' 
                    is not installed!")
                }else{
                    library(parallel)
                    if(is.null(nrCores)){
                        nrCores <- detectCores()
                    }
                }
            }else{
                useMulticore <- FALSE
                warning("Multicore is not supported on non-Unix platforms! 
                Setting useMulticore=FALSE")
            }
        }
        
        if(!(length(unique(width(regionsSeq))) == 1)){
            stop("All sequences in the input DNAStringSet must have the 
            same length!")
        }
        if(!(length(seqOrder) == length(regionsSeq))){
            stop("The number of elements in 'seqOrder' must match the number 
            of input sequences in 'regionsSeq'!")
        }
        if(length(patterns) == 0){
            stop("At least one pattern needs to be specified!")
        }
        
        if(length(flankUp) == 0){
            flankUp <- round(width(regionsSeq)[1]/2)
        }
        if(length(flankDown) == 0){
            flankDown <- width(regionsSeq)[1] - flankUp
        }
        if(!(color %in% c("green","cyan","blue","purple","pink","red","orange",
        "brown","gray"))){
            stop("Specified color scale not supported! Please choose one of the 
following color scales: 'green','cyan','blue','purple','pink','red','orange',
'brown','gray', and refer to the vignette for their appearance.")
        }
        
        message("\nGetting oligonucleotide occurrence matrix...")
        patterns.occurence.melted.list <- getPatternOccurrenceList(regionsSeq =
        regionsSeq, patterns = patterns, seqOrder = seqOrder,
        useMulticore = useMulticore, nrCores = nrCores)
        
        a <- .pattern.smoothscatter(melted = patterns.occurence.melted.list,
        orig = regionsSeq, patterns = patterns, flankUp = flankUp,
        flankDown = flankDown, bw = bandWidth, nbin = nBin, color = color,
        transf = transf, xTicks = xTicks, xTicksAt = xTicksAt, xLabel = xLabel,
        yTicks = yTicks, yTicksAt = yTicksAt, yLabel = yLabel,
        cex.axis = cexAxis, plot.scale = plotScale, scale.length = scaleLength,
        scale.width = scaleWidth, add.label = addPatternLabel,
        cex.label = cexLabel, label.col = labelCol,
        addReferenceLine = addReferenceLine, plotColorLegend = plotColorLegend,
        out = outFile, plot.width = plotWidth, plot.height = plotHeight,
        useMulticore = useMulticore, nrCores = nrCores)
    
    }
)


setGeneric(
name="plotPatternOccurrenceAverage",
def=function(regionsSeq, patterns, flankUp = NULL, flankDown = NULL,
smoothingWindow = 1, color = rainbow(length(patterns)),
xLabel = "Distance to reference point (bp)", yLabel = "Relative frequency",
cexAxis = 1, addReferenceLine = TRUE, plotLegend = TRUE,
cexLegend = 1, useMulticore = FALSE, nrCores = NULL, add = FALSE, ...){
    standardGeneric("plotPatternOccurrenceAverage")
}
)

setMethod("plotPatternOccurrenceAverage",
signature(regionsSeq = "DNAStringSet"),
function(regionsSeq, patterns, flankUp = NULL, flankDown = NULL,
    smoothingWindow = 1, color = rainbow(length(patterns)),
    xLabel = "Distance to reference point (bp)", yLabel = "Relative frequency",
    cexAxis = 1, addReferenceLine = TRUE, plotLegend = TRUE,
    cexLegend = 1, useMulticore = FALSE, nrCores = NULL, add = FALSE, ...){
        
        if(!(length(unique(width(regionsSeq))) == 1)){
            stop("All sequences in the input DNAStringSet must have the
            same length!")
        }
        if(length(patterns) == 0){
            stop("At least one pattern needs to be specified!")
        }
        
        if(length(flankUp) == 0){
            flankUp <- round(width(regionsSeq)[1]/2)
        }
        if(length(flankDown) == 0){
            flankDown <- width(regionsSeq)[1] - flankUp
        }
        
        message("\nGetting oligonucleotide occurrence matrix...")
        patterns.occurence.melted.list <- getPatternOccurrenceList(regionsSeq =
        regionsSeq, patterns = patterns, seqOrder = c(1:length(regionsSeq)),
        useMulticore = useMulticore, nrCores = nrCores)
        pattern.widths <- width(patterns)
        names(pattern.widths) <- patterns
        
        .plot.windowed.average(occurence.melted.list =
        patterns.occurence.melted.list, nr.seq = length(regionsSeq),
        pattern.widths = pattern.widths, flankUp = flankUp, flankDown =
        flankDown, smoothingWindow = smoothingWindow, color = color, xLabel =
        xLabel, yLabel = yLabel, cexAxis = cexAxis, addReferenceLine =
        addReferenceLine, plotLegend = plotLegend, cexLegend = cexLegend,
        add = add, ...)

    }
)


setGeneric(
name="plotMotifDensityMap",
def=function(regionsSeq, motifPWM, minScore = "80%",
    seqOrder = c(1:length(regionsSeq)), flankUp = NULL, flankDown = NULL,
    nBin = NULL, bandWidth = NULL, color = "blue", transf = NULL, xTicks = NULL,
    xTicksAt = NULL, xLabel = "", yTicks = NULL, yTicksAt = NULL, yLabel = "",
    cexAxis = 8, plotScale = TRUE, scaleLength = NULL, scaleWidth = 15,
    addReferenceLine = TRUE, plotColorLegend = TRUE, outFile = "DensityMap",
    plotWidth = 2000, plotHeight = 2000){
        standardGeneric("plotMotifDensityMap")
    }
)

setMethod("plotMotifDensityMap",
signature(regionsSeq = "DNAStringSet", motifPWM = "matrix"),
function(regionsSeq, motifPWM, minScore = "80%",
    seqOrder = c(1:length(regionsSeq)), flankUp = NULL, flankDown = NULL,
    nBin = NULL, bandWidth = NULL, color = "blue", transf = NULL, xTicks = NULL,
    xTicksAt = NULL, xLabel = "", yTicks = NULL, yTicksAt = NULL, yLabel = "",
    cexAxis = 8, plotScale = TRUE, scaleLength = NULL, scaleWidth = 15,
    addReferenceLine = TRUE, plotColorLegend = TRUE, outFile = "DensityMap",
    plotWidth = 2000, plotHeight = 2000){
        
        if(!(length(unique(width(regionsSeq))) == 1)){
            stop("All sequences in the input DNAStringSet must have the 
            same length!")
        }
        if(!(length(seqOrder) == length(regionsSeq))){
            stop("The number of elements in 'seqOrder' must match the number 
            of input sequences in 'regionsSeq'!")
        }
        
        if(length(flankUp) == 0){
            flankUp <- round(width(regionsSeq)[1]/2)
        }
        if(length(flankDown) == 0){
            flankDown <- width(regionsSeq)[1] - flankUp
        }
        
        if(!(color %in% c("green","cyan","blue","purple","pink","red","orange",
        "brown","gray"))){
            stop("Specified color scale not supported! Please choose one of the 
following color scales: 'green','cyan','blue','purple','pink','red','orange','
brown','gray', and refer to the vignette for their appearance.")
        }

        message("\nGetting motif occurrence matrix...")
        motif.occurence.melted <- motifScanHits(regionsSeq = regionsSeq,
        motifPWM = motifPWM, minScore = minScore, seqOrder = seqOrder)
        
        if(nrow(motif.occurence.melted) == 0){
            cat("\nNo motif hits above specified score threshold found!\n
Try lowering the score threshold.\nExiting without making a density plot.\n\n")
        }else{
        
        motif.occurence.melted.list <- list(motif = motif.occurence.melted)
        
        a <- .pattern.smoothscatter(melted = motif.occurence.melted.list,
        orig = regionsSeq, patterns = "motif", flankUp = flankUp,
        flankDown = flankDown, bw = bandWidth, nbin = nBin, color = color,
        transf = transf, xTicks = xTicks, xTicksAt = xTicksAt, xLabel = xLabel,
        yTicks=yTicks, yTicksAt=yTicksAt, yLabel = yLabel, cex.axis = cexAxis,
        plot.scale = plotScale, scale.length = scaleLength,
        scale.width = scaleWidth, add.label = FALSE, cex.label = cexAxis,
        addReferenceLine = addReferenceLine, plotColorLegend = plotColorLegend,
        out = outFile, plot.width = plotWidth, plot.height = plotHeight)
        }
    }
)


setGeneric(
name="plotMotifOccurrenceAverage",
def=function(regionsSeq, motifPWM, minScore = "80%", flankUp = NULL,
flankDown = NULL, smoothingWindow = 1, color = "black",
xLabel = "Distance to reference point (bp)", yLabel = "Relative frequency",
cexAxis = 1, addReferenceLine = TRUE, plotLegend = FALSE, cexLegend = 1,
add = FALSE, ...){
    standardGeneric("plotMotifOccurrenceAverage")
}
)

setMethod("plotMotifOccurrenceAverage",
signature(regionsSeq = "DNAStringSet", motifPWM = "matrix"),
function(regionsSeq, motifPWM, minScore = "80%", flankUp = NULL,
flankDown = NULL, smoothingWindow = 1, color = "black",
xLabel = "Distance to reference point (bp)", yLabel = "Relative frequency",
cexAxis = 1, addReferenceLine = TRUE, plotLegend = FALSE, cexLegend = 1,
add = FALSE, ...){
    
    if(!(length(unique(width(regionsSeq))) == 1)){
        stop("All sequences in the input DNAStringSet must have the
        same length!")
    }
    if(length(flankUp) == 0){
        flankUp <- round(width(regionsSeq)[1]/2)
    }
    if(length(flankDown) == 0){
        flankDown <- width(regionsSeq)[1] - flankUp
    }
    
    message("\nGetting motif occurrence matrix...")
    motif.occurence.melted <- motifScanHits(regionsSeq = regionsSeq, motifPWM =
        motifPWM, minScore = minScore, seqOrder = c(1:length(regionsSeq)))
    motif.occurence.melted.list <- list(motif = motif.occurence.melted)
    pattern.widths <- ncol(motifPWM)
    names(pattern.widths) <- "motif"
    
    .plot.windowed.average(occurence.melted.list = motif.occurence.melted.list,
    nr.seq = length(regionsSeq), pattern.widths = pattern.widths, flankUp =
    flankUp, flankDown = flankDown, smoothingWindow = smoothingWindow, color =
    color, xLabel = xLabel, yLabel = yLabel, cexAxis = cexAxis,
    addReferenceLine = addReferenceLine, plotLegend = plotLegend, cexLegend =
    cexLegend, add = add, ...)
    
}
)


setGeneric(
name="plotMotifScanScores",
def=function(regionsSeq, motifPWM, seqOrder = c(1:length(regionsSeq)),
    flankUp = NULL, flankDown = NULL, xTicks = NULL, xTicksAt = NULL,
    xLabel = "", yTicks = NULL, yTicksAt = NULL, yLabel = "", cexAxis = 8,
    plotScale = TRUE, scaleLength = NULL, scaleWidth = 15,
    addReferenceLine = TRUE, plotColorLegend = TRUE,
    outFile = "MotifScanningScores.png", plotWidth = 2000, plotHeight = 2000){
        standardGeneric("plotMotifScanScores")
    }
)

setMethod("plotMotifScanScores",
signature(regionsSeq = "DNAStringSet", motifPWM = "matrix"),
function(regionsSeq, motifPWM, seqOrder = c(1:length(regionsSeq)),
    flankUp = NULL, flankDown = NULL, xTicks = NULL, xTicksAt = NULL,
    xLabel = "", yTicks = NULL, yTicksAt = NULL, yLabel = "", cexAxis = 8,
    plotScale = TRUE, scaleLength = NULL, scaleWidth = 15,
    addReferenceLine = TRUE, plotColorLegend = TRUE,
    outFile = "MotifScanningScores.png", plotWidth = 2000, plotHeight = 2000){
        
        if(!(length(unique(width(regionsSeq))) == 1)){
            stop("All sequences in the input DNAStringSet must have the 
            same length!")
        }
        if(!(length(seqOrder) == length(regionsSeq))){
            stop("The number of elements in 'seqOrder' must match the number 
            of input sequences in 'regionsSeq'!")
        }
        
        if(length(flankUp) == 0){
            flankUp <- round(width(regionsSeq)[1]/2)
        }
        if(length(flankDown) == 0){
            flankDown <- width(regionsSeq)[1] - flankUp
        }
        
        message("\nGetting motif scanning scores...")
        motif.scanning.scores <- motifScanScores(regionsSeq = regionsSeq,
        motifPWM = motifPWM, seqOrder = seqOrder)
        
        message("\nPlotting heatmap...")
        png(filename=outFile, width = plotWidth, height = plotHeight)
        
        breaks <- c(0,60,seq(62,100,2))
        colsFun <- colorRampPalette(c("darkblue", "blue", "cyan", "yellow",
        "red2"))
        cols <- colsFun(length(breaks) - 1)
        
        if(plotColorLegend){
            layout(mat = matrix(c(2,1), nrow = 1), widths = c(0.88, 0.12))
            par(mar = c(12, 0, 2, 12))
            plot(c(0,0), type = "n", xaxt = "n", yaxt = "n", xlab = "",
            ylab = "", xlim = c(0,1), ylim = c(0,1), xaxs = "i", yaxs = "i")
            box(lwd = 6)
            color.legend(xl=0, yb=0, xr=1, yt=1,
            legend=c("0 %", seq(20,100,20)), rect.col=rep(cols,
            times=round(diff(breaks)/min(diff(breaks)))),
            align="rb", gradient="y", cex=cexAxis)
        }
        
        .plot.motif.heatmap(motifScanningScores = motif.scanning.scores,
        flankUp = flankUp, flankDown = flankDown, cols = cols, breaks = breaks,
        xTicks = xTicks, xTicksAt = xTicksAt, xLabel = xLabel, yTicks = yTicks,
        yTicksAt = yTicksAt, yLabel = yLabel, cexAxis = cexAxis,
        plotScale = plotScale, scaleLength = scaleLength,
        scaleWidth = scaleWidth, addReferenceLine = addReferenceLine)
        
        dev.off()
    
    }
)


