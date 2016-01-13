setClass("PWM",
    representation=list(
        matrix="matrix",
        name="character",
        min.score="vector"),
    prototype=c(
        matrix=matrix(0, ncol=0, nrow=4),
        name="",
        min.score="80%"))

PWM = function(mat, name="", min.score="80%") {
    if (nrow(mat) != 4) stop("PWM@matrix must have 4 rows")
    new("PWM", mat, name, min.score)
}

