test_.get.pattern.occurence.melted <- function(){
    data(zebrafishPromoters)
    data(TBPpwm)
    checkEquals(seqPattern:::.get.pattern.occurence.melted(pattern = "TCTGA",
    seq = zebrafishPromoters[1:5], seqOrder = c(1,4,2,3,5)),
    data.frame(sequence = as.integer(c(1,2,3,4,4,4,4)),
    position = as.integer(c(565,987,311,120,372,592,930)), value = rep(1,7)))
}
