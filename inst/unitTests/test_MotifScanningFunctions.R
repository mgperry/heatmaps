test_.scan.sequence.with.pwm <- function(){
    data(zebrafishPromoters)
    data(TBPpwm)
    checkEquals(seqPattern:::.scan.sequence.with.pwm(pwm = TBPpwm,
    seq = zebrafishPromoters[[1]], minScore = "90%"), c(51, 53, 55, 204, 344, 905))
}

test_.get.scanning.score <- function(){
    data(zebrafishPromoters)
    data(TBPpwm)
    checkEqualsNumeric(seqPattern:::.get.scanning.score(pwm = TBPpwm,
    seq=zebrafishPromoters[[1]])[c(1:5)],
    c(-18.805607170,-23.526148870,-18.062235212,-23.998492469,-25.615083293),
    tolerance = 10^-9)
}

