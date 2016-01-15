#######################################
# Function for plotting average signal


plotMetaRegion <- function(occurence.melted.list, nr.seq,
pattern.widths, flankUp = NULL, flankDown = NULL, smoothingWindow = 3,
color = rainbow(length(occurence.melted.list)), xLabel =
"Distance to reference point (bp)", yLabel = "Relative frequency", cexAxis = 2,
addReferenceLine = TRUE, plotLegend = TRUE, cexLegend = 2, add = FALSE, ...){

    message("\nPlotting average signal...\n")
    patterns.avg.signal<-lapply(c(1:length(occurence.melted.list)),
    function(i) {
        x <- occurence.melted.list[[i]]
        a.s <- table(x$position)/nr.seq
        avg.signal <- rep(0, times =
        flankUp+flankDown-pattern.widths[i]+1)
        names(avg.signal) <- c(1:length(avg.signal))
        avg.signal[names(a.s)] <- a.s
        return(avg.signal)
    })

    starts <- lapply(patterns.avg.signal, function(x) {
        seq(1, length(x), by = smoothingWindow)
    })
    ends <- lapply(starts, function(x) {x + smoothingWindow - 1})
    ends <- lapply(1:length(ends), function(x) {
        ends[[x]][length(ends[[x]])] <- min(ends[[x]][length(ends[[x]])],
        length(patterns.avg.signal[[x]])); return(ends[[x]])
    })

    patterns.avg.signal.windowed <- lapply(c(1:length(patterns.avg.signal)),
    function(x){
        a.s <- sapply(seq(1:length(starts[[x]])), function(y) {
            mean(patterns.avg.signal[[x]][starts[[x]][y]:ends[[x]][y]])})
        x.coor <- (starts[[x]]+ends[[x]])/2 - flankUp - 1
        return(list(x.coor, a.s))
    })

    if(!(add)){
        plot(1, 1, xlim = c(-flankUp, flankDown), ylim = c(0,
        1.05*max(unlist(lapply(patterns.avg.signal.windowed, function(x) {
            max(x[[2]])})))), type = "n", xlab = xLabel, ylab = yLabel,
        cex.axis = cexAxis, cex.lab = cexAxis, ...)
    }

    a.s <- lapply(c(1:length(patterns.avg.signal.windowed)), function(i){
        lines(x = patterns.avg.signal.windowed[[i]][[1]],
        y = patterns.avg.signal.windowed[[i]][[2]], type = "l",
        col = color[i], ...)
    })

    if(plotLegend){
        legend("top", legend = names(pattern.widths), bty = "n", horiz = TRUE,
        lwd = 1, col = color, cex = cexLegend)
    }
    if(addReferenceLine){
        abline(v = 0, lty = "dashed")
    }


}

