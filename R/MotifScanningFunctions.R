.scan.sequence.with.pwm <- function(pwm, seq, minScore){
    
    nc <- nchar(minScore)
    if (substr(minScore, nc, nc) == "%"){
        perc.threshold <- substr(minScore, 1L, nc-1L)
        min.score <- minScore(pwm)
        max.score <- maxScore(pwm)
        score.threshold = min.score + as.double(perc.threshold)/100 *
            (max.score-min.score)
    }else{
        score.threshold <- minScore
    }

    # make all PWM scores positive to avoid matching the Ns which are assigned score 0
    pwm.corr <- pwm + abs(min(pwm))
    score.threshold <- score.threshold + ncol(pwm) * abs(min(pwm))

    pwm.match <- matchPWM(pwm = pwm.corr, subject = seq, min.score = score.threshold)
    return(start(pwm.match))

}

.get.scanning.score <- function(pwm, seq){

# make all PWM scores positive to avoid matching the Ns which are assigned score 0
    pwm.corr <- pwm + abs(min(pwm))

    pwm.scores <- PWMscoreStartingAt(pwm = pwm.corr, subject = seq,
    starting.at = c(1:(length(seq) - ncol(pwm) + 1)))

    pwm.scores <- pwm.scores - ncol(pwm) * abs(min(pwm))
    return(pwm.scores)

}

