### R code from vignette source 'seqPattern.Rnw'

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: setup
###################################################
options(width = 80)
olocale=Sys.setlocale(locale="C")


###################################################
### code chunk number 3: seqPattern.Rnw:93-94
###################################################
library(seqPattern)


###################################################
### code chunk number 4: seqPattern.Rnw:106-109
###################################################
data(zebrafishPromoters)
zebrafishPromoters
head(zebrafishPromoters@elementMetadata)


###################################################
### code chunk number 5: seqPattern.Rnw:125-127
###################################################
data(zebrafishPromoters24h)
head(zebrafishPromoters24h)


###################################################
### code chunk number 6: seqPattern.Rnw:150-151
###################################################
library(BSgenome.Drerio.UCSC.danRer7)


###################################################
### code chunk number 7: seqPattern.Rnw:158-161
###################################################
data(zebrafishPromoters24h)
nrow(zebrafishPromoters24h)
head(zebrafishPromoters24h)


###################################################
### code chunk number 8: seqPattern.Rnw:176-188
###################################################
# create GRanges of dominant TSS position with associated metadata
zebrafishPromotersTSS <- GRanges(seqnames = zebrafishPromoters24h$chr,
            ranges = IRanges(start = zebrafishPromoters24h$dominantTSS, 
            end = zebrafishPromoters24h$dominantTSS), 
            strand = zebrafishPromoters24h$strand, 
            interquantileWidth = zebrafishPromoters24h$interquantileWidth, 
            seqlengths = seqlengths(Drerio))
# get regions flanking TSS - 400 bp upstream and 800 bp downstream
zebrafishPromotersTSSflank <- promoters(zebrafishPromotersTSS, upstream = 400, 
                                        downstream = 800)
# obtain genomic sequence of flanking regions
zebrafishPromotersTSSflankSeq <- getSeq(Drerio, zebrafishPromotersTSSflank)


###################################################
### code chunk number 9: seqPattern.Rnw:204-208 (eval = FALSE)
###################################################
## plotPatternDensityMap(regionsSeq = zebrafishPromotersTSSflankSeq, 
##             patterns = c("AA", "TA", "CG", "GC"), 
##             seqOrder = order(zebrafishPromotersTSSflank$interquantileWidth), 
##             flankUp = 400, flankDown = 800, color = "blue")


###################################################
### code chunk number 10: seqPattern.Rnw:264-268 (eval = FALSE)
###################################################
## plotPatternDensityMap(regionsSeq = zebrafishPromotersTSSflankSeq, 
##             patterns = c("AA", "TA", "CG", "GC"), 
##             seqOrder = order(zebrafishPromotersTSSflank$interquantileWidth),
##             flankUp = 400, flankDown = 800, useMulticore = TRUE, nrCores = 4)


###################################################
### code chunk number 11: seqPattern.Rnw:310-314 (eval = FALSE)
###################################################
## plotPatternDensityMap(regionsSeq = zebrafishPromotersTSSflankSeq, 
##             patterns = c("WW", "SS"), 
##             seqOrder = order(zebrafishPromotersTSSflank$interquantileWidth), 
##             flankUp = 400, flankDown = 800, color = "cyan", labelCol = "white")


###################################################
### code chunk number 12: seqPattern.Rnw:333-345 (eval = FALSE)
###################################################
## # make index of sharp and broad promoters
## sIdx <- zebrafishPromotersTSSflank$interquantileWidth <= 9
## bIdx <- zebrafishPromotersTSSflank$interquantileWidth > 9
## # plot average dinucleotide profile for sharp promoters
## par(mfrow = c(1,2), mar = c(4.5,4,1,1))
## plotPatternOccurrenceAverage(regionsSeq = zebrafishPromotersTSSflankSeq[sIdx], 
##             patterns = c("WW", "SS"), flankUp = 400, flankDown = 800, 
##             smoothingWindow = 3, color = c("red3", "blue3"), cex.axis = 0.9)
## # plot average dinucleotide profile for broad promoters
## plotPatternOccurrenceAverage(regionsSeq = zebrafishPromotersTSSflankSeq[bIdx], 
##             patterns = c("WW", "SS"), flankUp = 400, flankDown = 800, 
##             smoothingWindow = 3, color = c("red3", "blue3"), cex.axis = 0.9)


###################################################
### code chunk number 13: seqPattern.Rnw:370-375
###################################################
# get regions flanking dominant TSS - 200bp upstream and downstream
zebrafishPromotersTSSflank <- promoters(zebrafishPromotersTSS, upstream = 200, 
                                        downstream = 200)
# obtain genomic sequence of the flanking regions
zebrafishPromotersTSSflankSeq <- getSeq(Drerio, zebrafishPromotersTSSflank)


###################################################
### code chunk number 14: seqPattern.Rnw:377-389 (eval = FALSE)
###################################################
## # plot density of TATA-box consensus sequence in pink
## plotPatternDensityMap(regionsSeq = zebrafishPromotersTSSflankSeq, 
##             patterns = c("TATAWAWR", "YWTWTATA"), 
##             seqOrder = order(zebrafishPromotersTSSflank$interquantileWidth), 
##             flankUp = 200, flankDown = 200, nBin = c(400, 10000), 
##             bandWidth = c(1,6), color = "pink", addPatternLabel = FALSE)
##  # plot density of GC-box consensus sequence in purple
## plotPatternDensityMap(regionsSeq = zebrafishPromotersTSSflankSeq, 
##             patterns = c("RGGMGGR", "YCCKCCY"), 
##             seqOrder = order(zebrafishPromotersTSSflank$interquantileWidth), 
##             flankUp = 200, flankDown = 200, nBin = c(400, 10000), 
##             bandWidth = c(1,6), color = "purple", addPatternLabel = FALSE)


###################################################
### code chunk number 15: seqPattern.Rnw:429-431
###################################################
data(TBPpwm)
TBPpwm


###################################################
### code chunk number 16: seqPattern.Rnw:439-443 (eval = FALSE)
###################################################
## plotMotifDensityMap(regionsSeq = zebrafishPromotersTSSflankSeq, 
##             motifPWM = TBPpwm, minScore = "90%", 
##             seqOrder = order(zebrafishPromotersTSSflank$interquantileWidth),
##             flankUp = 200, flankDown = 200, color = "red")


###################################################
### code chunk number 17: seqPattern.Rnw:450-454 (eval = FALSE)
###################################################
## plotMotifScanScores(regionsSeq = zebrafishPromotersTSSflankSeq, 
##             motifPWM = TBPpwm,
##             seqOrder = order(zebrafishPromotersTSSflank$interquantileWidth),
##             flankUp = 200, flankDown = 200)


###################################################
### code chunk number 18: seqPattern.Rnw:481-494 (eval = FALSE)
###################################################
## # make index of sharp and broad promoters
## sIdx <- zebrafishPromotersTSSflank$interquantileWidth <= 9
## bIdx <- zebrafishPromotersTSSflank$interquantileWidth > 9
## # plot average motif occurrence profile for sharp promoters
## plotMotifOccurrenceAverage(regionsSeq = zebrafishPromotersTSSflankSeq[sIdx], 
##             motifPWM = TBPpwm, minScore = "90%", flankUp = 200, flankDown = 200, 
##             smoothingWindow = 3, color = c("red3"), cex.axis = 0.9)
## # add average motif occurrence profile for broad promoters to the existing plot
## plotMotifOccurrenceAverage(regionsSeq = zebrafishPromotersTSSflankSeq[bIdx], 
##             motifPWM = TBPpwm, minScore = "90%", flankUp = 200, flankDown = 200, 
##             smoothingWindow = 3, color = c("blue3"), add = TRUE)
## legend("topright", legend = c("sharp", "broad"), col = c("red3", "blue3"), 
## bty = "n", lwd = 1)


###################################################
### code chunk number 19: seqPattern.Rnw:519-524
###################################################
motifOccurrence <- motifScanHits(regionsSeq = 
            zebrafishPromotersTSSflankSeq[1:50], 
            motifPWM = TBPpwm, minScore = "90%", seqOrder = 
            order(zebrafishPromotersTSSflank$interquantileWidth[1:50]))
head(motifOccurrence)


###################################################
### code chunk number 20: seqPattern.Rnw:538-543
###################################################
scanScores <- motifScanScores(regionsSeq = zebrafishPromotersTSSflankSeq[1:50], 
            motifPWM = TBPpwm, seqOrder = 
            order(zebrafishPromotersTSSflank$interquantileWidth[1:50]))
dim(scanScores)
scanScores[1:6,1:6]


###################################################
### code chunk number 21: seqPattern.Rnw:552-553
###################################################
sessionInfo()


