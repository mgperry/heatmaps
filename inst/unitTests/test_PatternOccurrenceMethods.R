test_getPatternOccurrenceList <- function(){
    data(zebrafishPromoters)
    data(TBPpwm)
    checkEquals(getPatternOccurrenceList(regionsSeq = zebrafishPromoters[1:5],
    patterns = c("TCTGWA", "ASWDGTW"), seqOrder = c(1,4,2,3,5)),
    list(TCTGWA = data.frame(sequence = as.integer(c(2,4,4)),
    position = as.integer(c(61,120,930)), value = rep(1,3)),
    ASWDGTW = data.frame(sequence = as.integer(c(1,1,1,2,3,3,3,5)),
    position = as.integer(c(226,685,796,897,73,364,974,542)), value=rep(1,8))))
}
