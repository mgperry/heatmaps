.get.scanning.score <- function(pwm, seq){

# make all PWM scores positive to avoid matching the Ns which are assigned score 0
    pwm.corr <- pwm + abs(min(pwm))

    pwm.scores <- PWMscoreStartingAt(pwm = pwm.corr, subject = seq,
    starting.at = c(1:(length(seq) - ncol(pwm) + 1)))

    pwm.scores <- pwm.scores - ncol(pwm) * abs(min(pwm))
    return(pwm.scores)

}

