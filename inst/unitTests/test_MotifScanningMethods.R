test_motifScanHits <- function(){
    data(zebrafishPromoters)
    data(TBPpwm)
    checkEquals(motifScanHits(regionsSeq = zebrafishPromoters[1:5],
    motifPWM = TBPpwm, minScore = "95%", seqOrder = c(1,4,2,3,5)),
    data.frame(sequence = as.integer(c(1,1,1,2,4,4)),
    position = as.integer(c(51,53,55,271,761,815)), value = rep(1,6)))
}

test_motifScanScores <- function(){
    library(Biostrings)
    data(zebrafishPromoters)
    data(TBPpwm)
    checkEqualsNumeric(as.vector(motifScanScores(regionsSeq =
    subseq(zebrafishPromoters[1:5], start = 1, end = 10),
    motifPWM = TBPpwm, seqOrder = c(1,4,2,3,5))),
    c(47.16491,43.42071,50.86560,51.06228,53.98393,39.21015,56.48875,
    66.81754,73.50455,31.59049,48.41759,47.59350,49.04233,47.88635,65.04803),
    tolerance = 10^-5)
}

