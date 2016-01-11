\name{TBPpwm}
\docType{data}
\alias{TBPpwm}
\title{Position-weight matrix for TATA-box binding protein motif}
\description{
This is a \code{\link{matrix}} object representing a position-weight matrix 
(PWM) for TATA-box binding protein motif. The scores in the matrix correspond
to log2 ratio between probability of observing base b (rows) at position i
(columns) and background probability of base b. This PWM has been derived using 
position frequency matrix from Jaspar database. It is intended to be used as an 
input data in running examples from \code{\link{seqPattern}} package help pages.
}
\usage{data(TBPpwm)}
\value{A \code{\link{matrix}} object}
\format{A \code{\link{matrix}} object}
\keyword{datasets}